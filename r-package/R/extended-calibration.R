# Convexity-based null-angle integration. M_2(u) is convex in u. Upper
# chords and extrapolated-chord lower bounds give an explicit integration
# bracket, conditional on accurate floating-point evidence evaluations.

.evx_affine_area <- function(a, b, anchor_u, anchor_y, slope, floor = 0) {
  # Integral of min(1,max(floor,anchor_y+slope*(sqrt(2)*sin(v)-anchor_u))).
  if (floor >= 1) return(b - a)
  cuts <- c(a, b)
  if (slope != 0) {
    roots_u <- anchor_u + (c(floor, 1) - anchor_y) / slope
    inside <- is.finite(roots_u) & abs(roots_u) < sqrt(2)
    roots <- asin(roots_u[inside] / sqrt(2))
    cuts <- sort(unique(c(cuts, roots[roots > a & roots < b])))
  }
  value <- 0
  for (k in seq_len(length(cuts) - 1L)) {
    left <- cuts[k]
    right <- cuts[k + 1L]
    h <- right - left
    mid <- (left + right) / 2
    umid <- sqrt(2) * sin(mid)
    ymid <- anchor_y + slope * (umid - anchor_u)
    if (ymid <= floor) part <- floor * h
    else if (ymid >= 1) part <- h
    else {
      q <- h / 2
      sinc_minus_one <- if (abs(q) < 1e-3) {
        -q^2 / 6 + q^4 / 120 - q^6 / 5040
      } else sin(q) / q - 1
      average <- anchor_y + slope * ((umid - anchor_u) + umid * sinc_minus_one)
      part <- h * min(1, max(floor, average))
    }
    value <- value + part
  }
  value
}

.evx_cell_bounds <- function(cell, log_c) {
  angles <- cell$v
  parts <- cell$parts
  total <- .ev_logadd(parts[, 1], parts[, 2])
  u <- sqrt(2) * sin(angles)
  widths <- diff(angles)
  # Independent monotonicity of positive and negative effect contributions.
  floors <- exp(pmin(0, c(
    .ev_logadd(parts[1, 1], parts[2, 2]),
    .ev_logadd(parts[2, 1], parts[3, 2])) - log_c))
  upper_rectangle <- exp(pmin(0, c(max(total[1:2]), max(total[2:3])) - log_c))
  lower <- sum(widths * floors) / pi
  upper <- sum(widths * upper_rectangle) / pi
  safe <- all(is.finite(total)) && max(total - log_c) < 600 &&
    all(diff(u) > 1e-12) && all(widths > 1e-7)
  if (!safe) return(c(lower = lower, upper = upper))
  f <- exp(total - log_c)
  slopes <- diff(f) / diff(u)
  if (any(!is.finite(slopes))) return(c(lower = lower, upper = upper))
  # Upper chords on the two half cells.
  high <- (.evx_affine_area(angles[1], angles[2], u[1], f[1], slopes[1]) +
    .evx_affine_area(angles[2], angles[3], u[2], f[2], slopes[2])) / pi
  # The opposite-half chord, extrapolated from the midpoint, is a lower
  # supporting line by convexity. The component-cross floor is also valid.
  low <- (.evx_affine_area(angles[1], angles[2], u[2], f[2], slopes[2], floors[1]) +
    .evx_affine_area(angles[2], angles[3], u[2], f[2], slopes[1], floors[2])) / pi
  c(lower = max(lower, low), upper = min(upper, high))
}

.evx_general_calibrator <- function(prior, tolerance, max_intervals) {
  state <- new.env(parent = emptyenv())
  state$values <- new.env(parent = emptyenv())
  state$cells <- NULL
  at <- function(v) {
    key <- sprintf("%.17g", v)
    if (exists(key, state$values, inherits = FALSE)) {
      return(get(key, state$values, inherits = FALSE))
    }
    if (v == -pi / 2) {
      ans <- c(positive = -Inf, negative = Inf)
      if (prior$alternative == "greater") ans[2] <- -Inf
    } else if (v == pi / 2) {
      ans <- c(positive = Inf, negative = -Inf)
      if (prior$alternative == "less") ans[1] <- -Inf
    } else {
      ans <- .evx_components(sqrt(2) * sin(v), 2, 2 * cos(v)^2, prior)
    }
    assign(key, ans, state$values)
    ans
  }
  cell <- function(a, b) {
    v <- c(a, (a + b) / 2, b)
    list(v = v, parts = t(vapply(v, at, c(positive = 0, negative = 0))))
  }
  initialize <- function() {
    if (is.null(state$cells)) {
      edges <- seq(-pi / 2, pi / 2, length.out = 9L)
      state$cells <- lapply(seq_len(8L), function(k) cell(edges[k], edges[k + 1L]))
    }
  }
  evaluate_one <- function(lc) {
    if (lc == -Inf) return(c(p_eville = 1, lower = 1, gap = 0, intervals = 0))
    if (lc == Inf) return(c(p_eville = 0, lower = 0, gap = 0, intervals = 0))
    initialize()
    cells <- state$cells
    bounds <- t(vapply(cells, .evx_cell_bounds, c(lower = 0, upper = 0), log_c = lc))
    ville <- exp(min(0, -lc))
    repeat {
      low <- sum(bounds[, 1])
      high <- min(1, ville, sum(bounds[, 2]))
      if (low > high + 5e-9) {
        stop("Calibration bounds crossed; inspect evidence quadrature.", call. = FALSE)
      }
      low <- min(low, high) # Remove roundoff at equal bounds only.
      gap <- high - low
      if (gap <= tolerance || length(cells) >= max_intervals) break
      j <- which.max(bounds[, 2] - bounds[, 1])
      old <- cells[[j]]
      if (old$v[2] == old$v[1] || old$v[2] == old$v[3]) break
      left <- cell(old$v[1], old$v[2])
      right <- cell(old$v[2], old$v[3])
      cells[[j]] <- left
      cells[[length(cells) + 1L]] <- right
      bounds[j, ] <- .evx_cell_bounds(left, lc)
      bounds <- rbind(bounds, .evx_cell_bounds(right, lc))
    }
    state$cells <- cells
    if (gap > tolerance) {
      warning(sprintf(paste0("Calibration gap %.3g exceeds tolerance %.3g; ",
        "returning the conservative upper value."), gap, tolerance), call. = FALSE)
    }
    c(p_eville = high, lower = low, gap = gap, intervals = length(cells))
  }
  function(log_evidence) {
    unique_values <- sort(unique(log_evidence))
    computed <- t(vapply(unique_values, evaluate_one,
                         c(p_eville = 0, lower = 0, gap = 0, intervals = 0)))
    as.data.frame(computed[match(log_evidence, unique_values), , drop = FALSE])
  }
}

#' Construct the null calibration for a fixed prior and alternative
#'
#' @param prior A fixed object from ev_prior().
#' @param tolerance Requested absolute width of the null-angle integration
#'   bracket, default 1e-4. This does not bound floating-point evidence error.
#' @param max_intervals Maximum number of adaptive angle intervals.
#' @return A reusable ev_calibration object. Evaluate it with ev_pvalues().
#' @details Monitoring starts at n=2. General proper priors use convexity of
#'   the initial mixture as a function of the invariant statistic to bracket
#'   the clipped null expectation. The returned upper value is conservative
#'   for exact evidence evaluations. Analytic formulas are used for centered
#'   two-sided Gaussian priors and full/half Lebesgue measure. The proper-prior
#'   upper bound is also capped by the ordinary Ville bound. No fixed cutoff
#'   in a simulated null sample is used.
#' @export
ev_calibration <- function(prior = ev_prior(), tolerance = 1e-4,
                           max_intervals = 4096L) {
  .evx_check_prior(prior)
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0 || tolerance >= 1 ||
      length(max_intervals) != 1L || !is.finite(max_intervals) ||
      max_intervals < 8 || max_intervals != as.integer(max_intervals)) {
    stop("Require 0 < tolerance < 1 and integer max_intervals >= 8.", call. = FALSE)
  }
  analytic <- !prior$proper ||
    (prior$family == "gaussian" && prior$location == 0 &&
       prior$alternative == "two.sided")
  if (analytic) {
    evaluate <- function(log_evidence) {
      if (prior$proper) p <- .ev_g_gaussian(log_evidence, prior$scale)
      else {
        log_a <- 0.5 * log(pi)
        if (prior$alternative != "two.sided") log_a <- log_a - log(4)
        z <- exp(pmin(0, (log_a - log_evidence) / 2))
        p <- (2 / pi) * (asin(z) + z * sqrt(pmax(0, 1 - z^2)))
      }
      data.frame(p_eville = p, lower = p, gap = 0, intervals = 0)
    }
  } else evaluate <- .evx_general_calibrator(prior, tolerance, max_intervals)
  structure(list(prior = prior, n0 = 2L, tolerance = tolerance,
    max_intervals = max_intervals,
    method = if (analytic) "analytic" else "adaptive convex bounds",
    evaluate = evaluate), class = "ev_calibration")
}

#' Calibrate log evidence with its matching null calibration
#'
#' @param log_evidence Numeric vector of natural log evidence; infinities are
#'   allowed. Use cumulative maxima for running p-values.
#' @param calibration Object from ev_calibration().
#' @return A data frame containing p_eville (upper value), lower, gap, and
#'   intervals (number of adaptive angle intervals).
#' @export
ev_pvalues <- function(log_evidence, calibration) {
  if (!inherits(calibration, "ev_calibration") ||
      !is.numeric(log_evidence) || anyNA(log_evidence) || !length(log_evidence)) {
    stop("Supply numeric log evidence without NA and an ev_calibration object.",
         call. = FALSE)
  }
  result <- calibration$evaluate(log_evidence)
  rownames(result) <- NULL
  result
}

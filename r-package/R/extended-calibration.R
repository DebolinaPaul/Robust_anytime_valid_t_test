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
#' @param tolerance Requested absolute width of the numerical calibration
#'   bracket, default 1e-5 as in the Python paper simulation. For Monte Carlo
#'   this brackets the empirical mean, not the population expectation.
#' @param max_intervals Maximum number of adaptive intervals.
#' @param method "monte_carlo" (default) or "deterministic" for proper priors.
#'   Lebesgue measure always uses its separate analytic calibration.
#' @param null_t Optional finite vector of at least two independently generated
#'   null t statistics with one degree of freedom. It must be independent of
#'   target observations and held fixed during monitoring. The default uses
#'   all 100000 draws from [ev_null_sample()], matching the Python notebook.
#' @return A reusable ev_calibration object. Evaluate it with [ev_pvalues()].
#' @details Monitoring starts at n=2. Monte Carlo calibration implements the
#'   same empirical clipped expectation and convex numerical brackets as the
#'   Python notebook. Its numerical upper estimate is capped by the ordinary
#'   Ville bound. This cap does not remove Monte Carlo error or establish an
#'   exact anytime-valid guarantee for the finite-sample approximation.
#'   The Monte Carlo standard error is descriptive, not a simultaneous bound.
#'
#'   Deterministic calibration brackets the population null-angle integral,
#'   with an upper value conservative for exact evidence evaluations. Centered
#'   two-sided Gaussian and full/half Lebesgue formulas are analytic in that
#'   mode. Numerical brackets do not certify every floating-point operation.
#' @seealso [ev_null_sample()], [ev_ttest()], [ev_pvalues()]
#' @export
ev_calibration <- function(prior = ev_prior(), tolerance = 1e-5,
                           max_intervals = 4096L,
                           method = c("monte_carlo", "deterministic"),
                           null_t = NULL) {
  .evx_check_prior(prior)
  method <- match.arg(method)
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0 || tolerance >= 1 ||
      !is.numeric(max_intervals) || length(max_intervals) != 1L ||
      !is.finite(max_intervals) || max_intervals < 8 ||
      max_intervals > .Machine$integer.max ||
      max_intervals != as.integer(max_intervals)) {
    stop("Require 0 < tolerance < 1 and integer max_intervals >= 8.", call. = FALSE)
  }
  if (!is.null(null_t) && (method != "monte_carlo" || !prior$proper)) {
    stop("`null_t` is only used by Monte Carlo calibration of a proper prior.",
         call. = FALSE)
  }
  if (prior$proper && method == "monte_carlo") {
    is_reference <- is.null(null_t)
    if (is_reference) null_t <- ev_null_sample()
    mc <- .evx_mc_calibrator(prior, null_t, tolerance, max_intervals)
    return(structure(list(prior = prior, n0 = 2L, tolerance = tolerance,
      max_intervals = max_intervals, method = "monte_carlo",
      B = mc$B, null_t = as.numeric(null_t),
      null_sample_source = if (is_reference) "Python notebook: seed 20260510" else "supplied null_t",
      approximation = TRUE, evaluate = mc$evaluate), class = "ev_calibration"))
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
    method = if (analytic) "analytic" else "deterministic",
    B = 0L, approximation = FALSE,
    evaluate = evaluate), class = "ev_calibration")
}

#' Calibrate log evidence with its matching null calibration
#'
#' @param log_evidence Numeric vector of natural log evidence; infinities are
#'   allowed. Use cumulative maxima for running p-values.
#' @param calibration Object from ev_calibration().
#' @param relative_tolerance Optional relative numerical gap for Monte Carlo
#'   evaluation. NULL uses only the calibrator's absolute tolerance.
#' @return A data frame with p_eville, lower, gap and intervals. In Monte Carlo
#'   mode, p_eville is the empirical numerical upper value capped by Ville;
#'   lower is capped in the same way. Additional columns mc_estimate,
#'   mc_lower, mc_upper, mc_gap, mcse and B describe the uncapped empirical
#'   calibration and sample size. These brackets do not bound Monte Carlo
#'   error relative to the population calibration.
#' @export
ev_pvalues <- function(log_evidence, calibration, relative_tolerance = NULL) {
  if (!inherits(calibration, "ev_calibration") ||
      !is.numeric(log_evidence) || anyNA(log_evidence) || !length(log_evidence)) {
    stop("Supply numeric log evidence without NA and an ev_calibration object.",
         call. = FALSE)
  }
  if (!is.null(relative_tolerance) &&
      (!is.numeric(relative_tolerance) || length(relative_tolerance) != 1L ||
       !is.finite(relative_tolerance) || relative_tolerance <= 0 ||
       relative_tolerance >= 1)) {
    stop("`relative_tolerance` must be NULL or a number between zero and one.",
         call. = FALSE)
  }
  result <- if (identical(calibration$method, "monte_carlo")) {
    calibration$evaluate(log_evidence, relative_tolerance = relative_tolerance)
  } else {
    calibration$evaluate(log_evidence)
  }
  rownames(result) <- NULL
  result
}

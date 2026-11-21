# Monte Carlo calibration matched to MCCalibration in Code_simulation.ipynb.
# Convexity brackets the empirical mean for the fixed null sample. It does
# not bound Monte Carlo error relative to the population calibration.

.evx_mc_cumsum <- function(x) {
  # Compensated prefix sums reduce cancellation on R platforms where
  # long double has the same precision as double.
  answer <- numeric(length(x) + 1L)
  total <- 0
  correction <- 0
  for (i in seq_along(x)) {
    adjusted <- x[i] - correction
    next_total <- total + adjusted
    correction <- (next_total - total) - adjusted
    total <- next_total
    answer[i + 1L] <- total
  }
  answer
}

.evx_mc_calibrator <- function(prior, null_t, tolerance = 1e-6,
                               max_intervals = 4096L,
                               log_evidence_fn = NULL, log_mean = 0) {
  if (!is.numeric(null_t) || !is.null(dim(null_t)) || length(null_t) < 2L ||
      any(!is.finite(null_t))) {
    stop("Supply at least two finite numeric null t statistics.", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0 || tolerance >= 1 ||
      !is.numeric(max_intervals) || length(max_intervals) != 1L ||
      !is.finite(max_intervals) || max_intervals < 8 ||
      max_intervals != floor(max_intervals) ||
      !is.numeric(log_mean) || length(log_mean) != 1L || !is.finite(log_mean)) {
    stop("Require 0 < tolerance < 1, integer max_intervals >= 8 and finite log_mean.",
         call. = FALSE)
  }
  if (!is.null(log_evidence_fn) && !is.function(log_evidence_fn)) {
    stop("The log-evidence callback must be a function.", call. = FALSE)
  }
  t <- sort(as.numeric(null_t))
  B <- length(t)
  h <- pmax(1, abs(t))
  u <- sqrt(2) * (t / h) / sqrt((1 / h)^2 + (t / h)^2)
  prefix_u <- .evx_mc_cumsum(u)
  prefix_u2 <- .evx_mc_cumsum(u^2)
  state <- new.env(parent = emptyenv())
  state$values <- new.env(parent = emptyenv())
  state$cells <- NULL
  state$tail_cache <- new.env(parent = emptyenv())
  state$thresholds <- new.env(parent = emptyenv())
  state$max_gap <- 0

  unit_exp <- function(x) exp(min(0, x))
  # Sample indices below are zero-based, matching the notebook. Intervals
  # include their left endpoint and exclude the right endpoint.
  at <- function(index) {
    key <- as.character(index)
    if (exists(key, state$values, inherits = FALSE)) {
      return(get(key, state$values, inherits = FALSE))
    }
    value_t <- t[index + 1L]
    if (is.null(log_evidence_fn)) {
      scale <- max(1, abs(value_t))
      log_d <- log(2) - 2 * log(scale) -
        log((1 / scale)^2 + (value_t / scale)^2)
      st <- list(n = 2L, t = value_t, d = exp(log_d))
      if (st$d <= 0) {
        stop("A sampled null statistic exceeds the evidence engine's numerical range.",
             call. = FALSE)
      }
      value <- .evx_log_mixture_dispatch(st, prior)
    } else value <- log_evidence_fn(value_t)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop("Nonfinite evidence at a sampled null statistic.", call. = FALSE)
    }
    assign(key, as.numeric(value), state$values)
    as.numeric(value)
  }
  cell <- function(start, stop) {
    indices <- c(start, floor((start + stop) / 2), stop)
    list(indices = indices, logs = vapply(indices, at, numeric(1)))
  }
  centered_moments <- function(start, stop, center) {
    count <- stop - start
    if (count <= 0) return(c(count = 0, first = 0, second = 0))
    if (count <= 256 || u[stop] - u[start + 1L] < 1e-3) {
      delta <- u[seq.int(start + 1L, stop)] - center
      return(c(count = count, first = sum(delta), second = sum(delta^2)))
    }
    first <- prefix_u[stop + 1L] - prefix_u[start + 1L]
    second <- prefix_u2[stop + 1L] - prefix_u2[start + 1L]
    c(count = count, first = first - count * center,
      second = max(0, second - 2 * center * first + count * center^2))
  }
  constant_stats <- function(start, stop, lower, upper) {
    count <- stop - start
    midpoint <- 0.5 * (lower + upper)
    c(count * lower, count * upper, count * midpoint^2)
  }
  pair_stats <- function(start, stop, anchor, low_value, low_slope,
                         high_value, high_slope, threshold) {
    if (stop <= start) return(numeric(3))
    cuts <- c(start, stop)
    for (j in seq_len(2L)) {
      value <- c(low_value, high_value)[j]
      slope <- c(low_slope, high_slope)[j]
      if (slope == 0) next
      for (level in c(0, threshold)) {
        difference <- level - value
        if (abs(difference) > 3 * abs(slope)) next
        root <- anchor + difference / slope
        # Number of sample U values strictly smaller than root, equivalent
        # to numpy.searchsorted(u, root, side="left").
        k <- findInterval(root, u, left.open = TRUE)
        if (start < k && k < stop) cuts <- c(cuts, k)
      }
    }
    cuts <- sort(unique(cuts))
    total <- numeric(3)
    for (j in seq_len(length(cuts) - 1L)) {
      left <- cuts[j]
      right <- cuts[j + 1L]
      center <- u[floor((left + right - 1) / 2) + 1L]
      coefficients <- function(value, slope) {
        at_center <- value + slope * (center - anchor)
        if (at_center <= 0) return(c(0, 0))
        if (at_center >= threshold) return(c(1, 0))
        c(at_center / threshold, slope / threshold)
      }
      low <- coefficients(low_value, low_slope)
      high <- coefficients(high_value, high_slope)
      moments <- centered_moments(left, right, center)
      count <- moments[1L]
      first <- moments[2L]
      second <- moments[3L]
      lower <- count * low[1L] + low[2L] * first
      upper <- count * high[1L] + high[2L] * first
      mid <- (low + high) / 2
      moment2 <- count * mid[1L]^2 + 2 * mid[1L] * mid[2L] * first +
        mid[2L]^2 * second
      total <- total + pmin(count, pmax(0, c(lower, upper, moment2)))
    }
    unname(total)
  }
  cell_stats <- function(one_cell, log_c) {
    a <- one_cell$indices[1L]
    m <- one_cell$indices[2L]
    b <- one_cell$indices[3L]
    logs <- one_cell$logs
    values_u <- u[c(a, m, b) + 1L]
    if (values_u[1L] == values_u[3L] && t[a + 1L] != t[b + 1L]) {
      values <- vapply(seq.int(a, b - 1L), function(i) unit_exp(at(i) - log_c),
                       numeric(1))
      return(c(sum(values), sum(values), sum(values^2)))
    }
    if (b - a == 1 || values_u[1L] == values_u[3L]) {
      value <- unit_exp(logs[1L] - log_c)
      return(constant_stats(a, b, value, value))
    }
    if (values_u[1L] == values_u[2L] || values_u[2L] == values_u[3L]) {
      return(constant_stats(a, b, 0, unit_exp(max(logs) - log_c)))
    }
    shift <- max(logs, log_c)
    f <- exp(logs - shift)
    threshold <- exp(log_c - shift)
    slopes <- diff(f) / diff(values_u)
    if (any(!is.finite(slopes))) {
      stop("Nonfinite affine slope in Monte Carlo calibration.", call. = FALSE)
    }
    if (threshold < 1e-250) {
      left_floor <- min(f[2L], f[2L] + slopes[2L] * (values_u[1L] - values_u[2L]))
      right_floor <- min(f[2L], f[2L] + slopes[1L] * (values_u[3L] - values_u[2L]))
      left_lower <- as.numeric(left_floor > 0 && log(left_floor) + shift >= log_c)
      right_lower <- as.numeric(right_floor > 0 && log(right_floor) + shift >= log_c)
      return(constant_stats(a, m, left_lower, unit_exp(max(logs[1:2]) - log_c)) +
        constant_stats(m, b, right_lower, unit_exp(max(logs[2:3]) - log_c)))
    }
    pair_stats(a, m, values_u[2L], f[2L], slopes[2L],
               f[2L], slopes[1L], threshold) +
      pair_stats(m, b, values_u[2L], f[2L], slopes[1L],
                 f[2L], slopes[2L], threshold)
  }
  output <- function(estimate, lower, upper, gap, mcse, used, ville,
                      intervals) {
    capped_lower <- min(lower, ville)
    c(p_eville = used, lower = capped_lower,
      gap = max(0, used - capped_lower), intervals = intervals,
      mc_estimate = estimate, mc_lower = lower, mc_upper = upper,
      mc_gap = gap, mcse = mcse, B = B, ville_bound = ville)
  }
  evaluate_one <- function(log_c, relative_tolerance = NULL,
                            .absolute_tolerance = NULL) {
    if (!is.numeric(log_c) || length(log_c) != 1L || is.na(log_c)) {
      stop("Supply one log threshold without NA or NaN.", call. = FALSE)
    }
    if (!is.null(relative_tolerance) &&
        (!is.numeric(relative_tolerance) || length(relative_tolerance) != 1L ||
         !is.finite(relative_tolerance) || relative_tolerance <= 0 ||
         relative_tolerance >= 1)) {
      stop("relative_tolerance must lie between zero and one.", call. = FALSE)
    }
    absolute <- if (is.null(.absolute_tolerance)) tolerance else .absolute_tolerance
    if (!is.finite(log_c)) {
      value <- as.numeric(log_c < 0)
      return(output(value, value, value, 0, 0, value, value, 0))
    }
    endpoint_max <- max(at(0L), at(B - 1L))
    if (log_c > endpoint_max) {
      factor <- exp(endpoint_max - log_c)
      anchor_tolerance <- exp(min(0, log(absolute) + log_c - endpoint_max))
      key <- if (is.null(relative_tolerance)) "none" else sprintf("%.17g", relative_tolerance)
      cached <- if (exists(key, state$tail_cache, inherits = FALSE)) {
        get(key, state$tail_cache, inherits = FALSE)
      } else NULL
      if (is.null(cached) || cached[["mc_gap"]] > anchor_tolerance) {
        cached <- evaluate_one(endpoint_max, relative_tolerance,
                                .absolute_tolerance = anchor_tolerance)
        assign(key, cached, state$tail_cache)
      }
      lower <- cached[["mc_lower"]] * factor
      upper <- cached[["mc_upper"]] * factor
      ville <- unit_exp(log_mean - log_c)
      gap <- cached[["mc_gap"]] * factor
      if (is.null(.absolute_tolerance)) state$max_gap <- max(state$max_gap, gap)
      return(output(cached[["mc_estimate"]] * factor, lower, upper, gap,
        cached[["mcse"]] * factor, min(upper, ville), ville, cached[["intervals"]]))
    }
    if (is.null(state$cells)) {
      edges <- unique(floor(seq(0, B - 1L, length.out = 9L)))
      state$cells <- lapply(seq_len(length(edges) - 1L), function(i) {
        cell(edges[i], edges[i + 1L])
      })
    }
    cells <- state$cells
    summaries <- lapply(cells, cell_stats, log_c = log_c)
    last <- unit_exp(at(B - 1L) - log_c)
    total <- colSums(do.call(rbind, summaries)) + c(last, last, last^2)
    repeat {
      lower <- total[1L] / B
      upper <- total[2L] / B
      if (lower > upper + 1e-9) {
        stop("Empirical calibration bounds crossed.", call. = FALSE)
      }
      lower <- min(lower, upper)
      gap <- max(0, upper - lower)
      requested <- absolute
      if (!is.null(relative_tolerance) && upper > 0) {
        requested <- min(requested, relative_tolerance * upper)
      }
      if (gap <= requested) break
      if (length(cells) >= max_intervals) {
        stop(sprintf(paste0("Empirical calibration gap %.3g exceeds %.3g; ",
                            "increase max_intervals."), gap, requested), call. = FALSE)
      }
      j <- which.max(vapply(summaries, function(x) x[2L] - x[1L], numeric(1)))
      indices <- cells[[j]]$indices
      a <- indices[1L]
      m <- indices[2L]
      b <- indices[3L]
      if (m == a || m == b) {
        stop("Sample interval cannot be refined further.", call. = FALSE)
      }
      left <- cell(a, m)
      right <- cell(m, b)
      ls <- cell_stats(left, log_c)
      rs <- cell_stats(right, log_c)
      total <- total + ls + rs - summaries[[j]]
      cells[[j]] <- left
      summaries[[j]] <- ls
      cells[[length(cells) + 1L]] <- right
      summaries[[length(summaries) + 1L]] <- rs
    }
    state$cells <- cells
    estimate <- 0.5 * (lower + upper)
    second <- total[3L] / B
    mcse <- sqrt(max(0, second - estimate^2) / (B - 1))
    ville <- unit_exp(log_mean - log_c)
    if (is.null(.absolute_tolerance)) state$max_gap <- max(state$max_gap, gap)
    output(estimate, lower, upper, gap, mcse, min(upper, ville), ville, length(cells))
  }
  evaluate <- function(log_evidence, relative_tolerance = NULL) {
    if (!is.numeric(log_evidence) || anyNA(log_evidence) || !length(log_evidence)) {
      stop("Supply numeric log evidence without NA or NaN.", call. = FALSE)
    }
    values <- sort(unique(log_evidence))
    computed <- do.call(rbind, lapply(values, evaluate_one,
                                    relative_tolerance = relative_tolerance))
    answer <- as.data.frame(computed[match(log_evidence, values), , drop = FALSE])
    rownames(answer) <- NULL
    answer
  }
  log_upper <- function(log_c, relative_tolerance = NULL) {
    result <- evaluate_one(log_c, relative_tolerance)
    endpoint_max <- max(at(0L), at(B - 1L))
    if (is.finite(log_c) && log_c > endpoint_max) {
      key <- if (is.null(relative_tolerance)) "none" else sprintf("%.17g", relative_tolerance)
      anchor <- get(key, state$tail_cache, inherits = FALSE)
      return(min(log(anchor[["mc_upper"]]) + endpoint_max - log_c,
                 log_mean - log_c, 0))
    }
    log(result[["p_eville"]])
  }
  log_threshold <- function(alpha, log_tolerance = 2e-7) {
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
        alpha <= 0 || alpha >= 1 || !is.numeric(log_tolerance) ||
        length(log_tolerance) != 1L || !is.finite(log_tolerance) || log_tolerance <= 0) {
      stop("Require 0 < alpha < 1 and positive log_tolerance.", call. = FALSE)
    }
    key <- paste(sprintf("%.17g", alpha), sprintf("%.17g", log_tolerance), sep = ":")
    if (exists(key, state$thresholds, inherits = FALSE)) {
      return(get(key, state$thresholds, inherits = FALSE))
    }
    high <- log_mean - log(alpha)
    # This upper endpoint is valid because the empirical value is capped
    # by the ordinary Ville bound; add a floating-point rounding margin.
    high <- high + 4 * .Machine$double.eps * max(1, abs(high))
    low <- high - 2
    step <- 2
    bracketed <- FALSE
    for (i in seq_len(80L)) {
      if (evaluate_one(low)[["p_eville"]] > alpha) {
        bracketed <- TRUE
        break
      }
      high <- low
      step <- 2 * step
      low <- low - step
    }
    if (!bracketed) stop("Could not bracket empirical log threshold.", call. = FALSE)
    for (i in seq_len(80L)) {
      middle <- (low + high) / 2
      if (evaluate_one(middle)[["p_eville"]] <= alpha) high <- middle else low <- middle
      if (high - low <= log_tolerance) break
    }
    if (evaluate_one(high)[["p_eville"]] > alpha * (1 + 1e-12)) {
      stop("Returned empirical threshold failed its check.", call. = FALSE)
    }
    assign(key, high, state$thresholds)
    high
  }
  list(evaluate = evaluate, evaluate_one = evaluate_one,
       log_upper = log_upper, log_threshold = log_threshold,
       B = B, null_t_sorted = t,
       dkw_radius = function(delta = 0.01) {
         if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) ||
             delta <= 0 || delta >= 1) stop("Require 0 < delta < 1.", call. = FALSE)
         sqrt(log(2 / delta) / (2 * B))
       },
       diagnostics = function() list(max_empirical_gap = state$max_gap,
         evidence_evaluations = length(ls(state$values, all.names = TRUE))))
}

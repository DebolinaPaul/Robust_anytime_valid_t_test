# Full-line normal-scale mixtures matching Code_simulation.ipynb.
# Order refinement and omitted-tail bounds are numerical checks, not
# interval-arithmetic certificates for every floating-point operation.

.evx_invariant_u <- function(st) {
  h <- max(abs(st$t), sqrt(st$n - 1))
  sqrt(st$n) * (st$t / h) /
    sqrt((st$t / h)^2 + (sqrt(st$n - 1) / h)^2)
}

.evx_log_j_positive <- function(z, n) {
  previous <- 0.5 * log(2 * pi) + stats::pnorm(z, log.p = TRUE)
  current <- .ev_logadd(-0.5 * z^2, log(z) + previous)
  if (n > 2L) {
    for (k in seq_len(n - 2L)) {
      next_value <- .ev_logadd(log(z) + current, log(k) + previous)
      previous <- current
      current <- next_value
    }
  }
  current
}

.evx_log_i_negative <- function(a, n, order = 128L) {
  if (n == 2L) return(.evx_log_i2_negative(a))
  q <- .evx_gamma_quad(n, order)
  answer <- numeric(length(a))
  starts <- seq.int(1L, length(a), by = 2048L)
  for (start in starts) {
    ix <- seq.int(start, min(length(a), start + 2047L))
    aa <- a[ix]
    rate <- pmax(aa, sqrt(n))
    values <- outer(1 - aa / rate, q$x) -
      outer(0.5 / rate^2, q$x^2)
    values <- sweep(values, 2L, q$log_w, "+")
    maximum <- apply(values, 1L, max)
    answer[ix] <- lgamma(n) - n * log(rate) + maximum +
      log(rowSums(exp(values - maximum)))
  }
  answer
}

.evx_log_gaussian_logsd <- function(st, location, log_sd,
                                   radial_order = 128L) {
  n <- st$n
  u <- .evx_invariant_u(st)
  log_d <- log(st$d)
  log_a <- .ev_logadd(0, log(n) + 2 * log_sd)
  log_D <- .ev_logadd(0, log_d + 2 * log_sd)
  z <- location * u * exp(-0.5 * (log_a + log_D))
  prefactor <- 0.5 * (n - 1) * log_a - 0.5 * n * log_D
  constant <- (n / 2 - 1) * log(2) + lgamma(n / 2)
  if (location * u >= 0) {
    # Algebraically cancel the competing quadratic exponent terms.
    prefactor - 0.5 * location^2 * exp(log_d - log_D) +
      .evx_log_j_positive(z, n) - constant
  } else {
    prefactor - 0.5 * n * location^2 * exp(-log_a) +
      .evx_log_i_negative(-z, n, radial_order) - constant
  }
}

.evx_log_scale_weight <- function(x, family) {
  switch(family,
    cauchy = 0.5 * log(2 / pi) - x - 0.5 * exp(-2 * x),
    student_t5 = log(2) + 2.5 * log(2.5) - lgamma(2.5) -
      5 * x - 2.5 * exp(-2 * x),
    laplace = 2 * x - 0.5 * exp(2 * x),
    horseshoe = log(2 / pi) + x - .ev_logadd(0, 2 * x),
    stop("This family is not a normal-scale mixture.", call. = FALSE)
  )
}

.evx_scale_integral <- function(st, prior, lower, upper, order) {
  count <- max(1L, ceiling((upper - lower) / 0.5))
  edges <- seq(lower, upper, length.out = count + 1L)
  half <- diff(edges) / 2
  mid <- (edges[-1L] + edges[-length(edges)]) / 2
  q <- .evx_legendre(order)
  x <- as.vector(outer(half, q$x) + mid)
  weights <- as.vector(outer(log(half), q$log_w, "+"))
  log_density <- .evx_log_scale_weight(x, prior$family)
  live <- is.finite(log_density)
  values <- .evx_log_gaussian_logsd(st, prior$location,
    log(prior$scale) + x[live], radial_order = max(64L, 4L * order))
  answer <- .evx_logsum(values + log_density[live] + weights[live])
  if (!is.finite(answer)) {
    stop("The normal-scale mixture is outside numerical range.", call. = FALSE)
  }
  list(log_value = answer, nodes = sum(live))
}

.evx_scale_tail_logs <- function(st, prior, lower, upper) {
  n <- st$n
  log_d <- log(st$d)
  log_C <- 0.5 * n * (log(2 * n) - log_d)
  family <- prior$family
  lower_mass <- switch(family,
    horseshoe = min(0, log(2 / pi) + lower),
    laplace = min(0, 2 * lower - log(2)),
    cauchy = if (-lower > 350) -Inf else
      min(0, log(2) + stats::pnorm(-exp(-lower), log.p = TRUE)),
    student_t5 = if (lower >= 0) 0 else if (-2 * lower > 700) -Inf else
      2.5 * (1 - 2 * lower - exp(-2 * lower)),
    stop("Unknown normal-scale-mixture family.", call. = FALSE)
  )
  log_factor <- .ev_log_lebesgue(n, st$d) -
    0.5 * log(2 * pi) - log(prior$scale)
  upper_moment <- switch(family,
    horseshoe = -log(pi) - 2 * upper,
    cauchy = 0.5 * log(2 / pi) - log(2) - 2 * upper,
    student_t5 = log(2) + 2.5 * log(2.5) - lgamma(2.5) -
      log(6) - 6 * upper,
    laplace = if (upper > 350) -Inf else
      0.5 * log(2 * pi) + stats::pnorm(-exp(upper), log.p = TRUE)
  )
  c(lower = log_C + lower_mass, upper = log_factor + upper_moment)
}

.evx_full_mixture_details <- function(st, prior, rtol = 2e-9) {
  if (prior$family == "gaussian") {
    answer <- .evx_log_gaussian_logsd(st, prior$location, log(prior$scale), 128L)
    discrepancy <- 0
    if (prior$location * st$t < 0 && st$n > 2L) {
      fine <- .evx_log_gaussian_logsd(st, prior$location, log(prior$scale), 256L)
      discrepancy <- abs(fine - answer)
      if (discrepancy > 2e-10) {
        stop("Gaussian radial quadrature failed order refinement.", call. = FALSE)
      }
      answer <- fine
    }
    if (!is.finite(answer)) stop("Gaussian mixture exceeded numerical range.", call. = FALSE)
    return(list(log_value = answer, relative_refinement_change = discrepancy,
      relative_omitted_tail_bound = 0, log_scale_lower = NA_real_,
      log_scale_upper = NA_real_, order = 0L, nodes = 0L))
  }
  lower <- -6
  upper <- 6
  for (attempt in seq_len(12L)) {
    coarse <- .evx_scale_integral(st, prior, lower, upper, 16L)
    fine <- .evx_scale_integral(st, prior, lower, upper, 32L)
    discrepancy <- abs(expm1(coarse$log_value - fine$log_value))
    order <- 32L
    if (!is.finite(discrepancy) || discrepancy > rtol / 4) {
      refined <- .evx_scale_integral(st, prior, lower, upper, 64L)
      discrepancy <- abs(expm1(fine$log_value - refined$log_value))
      fine <- refined
      order <- 64L
      if (!is.finite(discrepancy) || discrepancy > rtol / 4) {
        stop("Normal-scale-mixture quadrature failed order refinement.", call. = FALSE)
      }
    }
    tails <- .evx_scale_tail_logs(st, prior, lower, upper)
    log_tail <- .ev_logadd(tails[1L], tails[2L])
    if (log_tail - fine$log_value <= log(rtol / 4)) {
      return(list(log_value = fine$log_value,
        relative_refinement_change = discrepancy,
        relative_omitted_tail_bound = exp(log_tail - fine$log_value),
        log_scale_lower = lower, log_scale_upper = upper,
        order = order, nodes = fine$nodes))
    }
    target <- fine$log_value + log(rtol / 16)
    log_C <- 0.5 * st$n * (log(2 * st$n) - log(st$d))
    if (tails[1L] > target) {
      lower <- switch(prior$family,
        horseshoe = min(lower - 2, target - log_C - log(2 / pi)),
        laplace = min(lower - 2, 0.5 * (target - log_C + log(2))),
        lower - 2)
    }
    if (tails[2L] > target) {
      slope <- if (prior$family == "student_t5") 6 else 2
      upper <- upper + max(2, (tails[2L] - target) / slope)
    }
    if (!is.finite(lower + upper) || upper - lower > 4000) {
      stop("Scale separation exceeds the supported numerical range.", call. = FALSE)
    }
  }
  stop("Normal-scale-mixture omitted-tail bounds did not converge.", call. = FALSE)
}

.evx_log_full_mixture_stats <- function(st, prior, rtol = 2e-9) {
  .evx_full_mixture_details(st, prior, rtol)$log_value
}

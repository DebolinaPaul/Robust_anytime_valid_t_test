# Exact density conventions used in the Python simulation notebook.
# The horseshoe is the normal / half-Cauchy scale mixture, not a logarithmic
# density approximation. All helpers are internal and use base R only.

.evx_scaled_e1_log <- function(log_a) {
  # log(exp(a) * E1(a)), accepting log(a) to retain tiny nonzero arguments.
  out <- numeric(length(log_a))
  zero <- log_a == -Inf
  out[zero] <- Inf
  small <- is.finite(log_a) & log_a <= 0
  if (any(small)) {
    la <- log_a[small]
    a <- exp(la)
    term <- -a
    series <- term
    for (k in seq.int(2L, 100L)) {
      term <- -term * a / k
      addition <- term / k
      series <- series + addition
      if (all(abs(addition) <= 2e-16 * pmax(1, abs(series)))) break
    }
    e1 <- -0.577215664901532860606512090082402431 - la - series
    if (any(e1 <= 0) || any(!is.finite(e1))) {
      stop("Horseshoe density series is outside numerical range.", call. = FALSE)
    }
    out[small] <- a + log(e1)
  }
  middle <- is.finite(log_a) & log_a > 0 & log_a < log(50)
  if (any(middle)) {
    # Continued fraction for the scaled exponential integral. Evaluating
    # E1(a) first and multiplying by exp(a) would overflow at large a.
    a <- exp(log_a[middle])
    b <- a + 1
    c <- rep(1e300, length(a))
    d <- 1 / b
    h <- d
    converged <- FALSE
    for (k in seq_len(1000L)) {
      ak <- -as.double(k)^2
      b <- b + 2
      denominator <- ak * d + b
      denominator[abs(denominator) < 1e-300] <- 1e-300
      d <- 1 / denominator
      c <- b + ak / c
      c[abs(c) < 1e-300] <- 1e-300
      multiplier <- c * d
      h <- h * multiplier
      if (all(abs(multiplier - 1) <= 2e-15)) {
        converged <- TRUE
        break
      }
    }
    if (!converged || any(h <= 0) || any(!is.finite(h))) {
      stop("Horseshoe density continued fraction did not converge.", call. = FALSE)
    }
    out[middle] <- log(h)
  }
  large <- log_a >= log(50)
  if (any(large)) {
    la <- log_a[large]
    inv <- exp(-la)
    term <- rep(1, length(la))
    series <- term
    for (k in seq_len(25L)) {
      term <- -term * k * inv
      series <- series + term
      if (all(abs(term) <= 2e-15 * abs(series))) break
    }
    if (any(series <= 0)) {
      stop("Horseshoe density tail expansion failed.", call. = FALSE)
    }
    out[large] <- -la + log(series)
  }
  if (anyNA(log_a)) stop("Horseshoe density received missing input.", call. = FALSE)
  out
}

.evx_prior_logpdf <- function(theta, family, location, scale) {
  if (family == "gaussian") {
    return(stats::dnorm(theta, mean = location, sd = scale, log = TRUE))
  }
  if (family == "cauchy") {
    return(stats::dcauchy(theta, location = location, scale = scale, log = TRUE))
  }
  if (family == "student_t5") {
    return(stats::dt((theta - location) / scale, df = 5, log = TRUE) - log(scale))
  }
  if (family == "laplace") {
    return(-log(2) - log(scale) - abs((theta - location) / scale))
  }
  if (family != "horseshoe") stop("Unknown proper prior family.", call. = FALSE)
  difference <- abs(theta - location)
  log_difference <- log(difference)
  overflow <- is.infinite(difference) & is.finite(theta) & is.finite(location)
  if (any(overflow)) {
    # Oppositely signed finite endpoints can overflow when subtracted.
    a <- log(abs(theta[overflow]))
    b <- log(abs(location))
    m <- pmax(a, b)
    log_difference[overflow] <- m + log1p(exp(-abs(a - b)))
  }
  log_a <- 2 * (log_difference - log(scale)) - log(2)
  .evx_scaled_e1_log(log_a) - log(scale) - log(pi) - 0.5 * log(2 * pi)
}

.evx_density_log_atan_exp <- function(x) {
  out <- x
  mid <- x >= -20 & x < 0
  out[mid] <- log(atan(exp(x[mid])))
  high <- x >= 0
  out[high] <- log(atan2(1, exp(-x[high])))
  out
}

.evx_horseshoe_log_sf_log_y <- function(log_y) {
  # For H = Z * lambda, Z ~ N(0,1), lambda ~ half-Cauchy(0,1),
  # P(H > y) = E[atan(abs(Z)/y)] / pi, y > 0. Scaling the integrand
  # by its value at abs(Z)=1 resolves the survival probability even when
  # y is enormous. Direct integration in H-space can miss this tail.
  if (log_y == -Inf) return(-log(2))
  if (!is.finite(log_y)) stop("Horseshoe tail argument is not finite.", call. = FALSE)
  anchor <- .evx_density_log_atan_exp(-log_y)
  result <- stats::integrate(function(z) {
    exp(-0.5 * z^2 + .evx_density_log_atan_exp(log(z) - log_y) - anchor)
  }, 0, Inf, rel.tol = 1e-11, abs.tol = 0, subdivisions = 1000L,
  stop.on.error = TRUE)
  if (!is.finite(result$value) || result$value <= 0 ||
      result$abs.error > 1e-9 * result$value) {
    stop("Horseshoe half-line probability failed its relative-error check.", call. = FALSE)
  }
  min(-log(2), 0.5 * log(2 / pi) - log(pi) + anchor + log(result$value))
}

.evx_prior_log_masses <- function(family, location, scale) {
  if (location == 0) return(c(positive = -log(2), negative = -log(2)))
  y <- abs(location) / scale
  if (family == "gaussian") {
    small <- stats::pnorm(y, lower.tail = FALSE, log.p = TRUE)
  } else if (family == "cauchy") {
    small <- stats::pcauchy(y, lower.tail = FALSE, log.p = TRUE)
  } else if (family == "student_t5") {
    small <- stats::pt(y, df = 5, lower.tail = FALSE, log.p = TRUE)
  } else if (family == "laplace") {
    small <- -log(2) - y
  } else if (family == "horseshoe") {
    small <- .evx_horseshoe_log_sf_log_y(log(abs(location)) - log(scale))
  } else stop("Unknown proper prior family.", call. = FALSE)
  large <- log(-expm1(small))
  if (location > 0) c(positive = large, negative = small)
  else c(positive = small, negative = large)
}

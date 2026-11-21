# Shifted Gaussian/Cauchy invariant mixtures. No noncentral-t density is used.
# Quadrature refinement is a numerical diagnostic, not interval arithmetic.

.evx_quad_cache <- new.env(parent = emptyenv())

.evx_logsum <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

.evx_legendre <- function(order) {
  key <- paste0("legendre_", order)
  if (exists(key, .evx_quad_cache, inherits = FALSE)) {
    return(get(key, .evx_quad_cache, inherits = FALSE))
  }
  k <- seq_len(order - 1L)
  off <- k / sqrt(4 * k^2 - 1)
  mat <- matrix(0, order, order)
  mat[cbind(k, k + 1L)] <- off
  mat[cbind(k + 1L, k)] <- off
  eig <- eigen(mat, symmetric = TRUE)
  ans <- list(x = eig$values, log_w = log(2 * eig$vectors[1L, ]^2))
  assign(key, ans, .evx_quad_cache)
  ans
}

.evx_gamma_quad <- function(n, order) {
  key <- paste0("gamma_", n, "_", order)
  if (exists(key, .evx_quad_cache, inherits = FALSE)) {
    return(get(key, .evx_quad_cache, inherits = FALSE))
  }
  k <- seq_len(order - 1L)
  mat <- diag(2 * seq_len(order) + n - 2)
  off <- sqrt(k * (k + n - 1))
  mat[cbind(k, k + 1L)] <- off
  mat[cbind(k + 1L, k)] <- off
  eig <- eigen(mat, symmetric = TRUE)
  ans <- list(x = eig$values, log_w = log(eig$vectors[1L, ]^2))
  assign(key, ans, .evx_quad_cache)
  ans
}

.evx_log_i2_negative <- function(a) {
  out <- numeric(length(a))
  small <- a < 10
  if (any(small)) {
    aa <- a[small]
    lr <- log(aa) + stats::pnorm(-aa, log.p = TRUE) +
      aa^2 / 2 + 0.5 * log(2 * pi)
    out[small] <- log(-expm1(lr))
  }
  if (any(!small)) {
    aa <- a[!small]
    inv <- 1 / aa^2
    term <- rep(1, length(aa))
    value <- term
    # Alternating expansion after integrating the Taylor expansion of
    # exp(-r^2/2). The next term bounds the integrated remainder.
    for (j in seq_len(40L)) {
      term <- -term * (2 * j + 1) * inv
      value <- value + term
      if (all(abs(term) <= 2e-16 * abs(value))) break
    }
    if (any(value <= 0) || any(abs(term) > 1e-12 * abs(value))) {
      stop("The negative radial tail failed its remainder check.", call. = FALSE)
    }
    out[!small] <- -2 * log(aa) + log(value)
  }
  out
}

.evx_log_lr <- function(theta, u, n, d, order = 64L) {
  if (any(!is.finite(theta)) || any(theta < 0)) {
    stop("Internal effects must be finite and nonnegative.", call. = FALSE)
  }
  if (u == 0) return(-n * theta^2 / 2)
  z <- u * theta
  chi_constant <- (n / 2 - 1) * log(2) + lgamma(n / 2)
  if (u > 0) {
    j1 <- 0.5 * log(2 * pi) + stats::pnorm(z, log.p = TRUE)
    j2 <- .ev_logadd(-z^2 / 2, log(z) + j1)
    if (n > 2L) {
      for (k in seq_len(n - 2L)) {
        next_j <- .ev_logadd(log(z) + j2, log(k) + j1)
        j1 <- j2
        j2 <- next_j
      }
    }
    ans <- -d * theta^2 / 2 + j2 - chi_constant
  } else if (n == 2L) {
    ans <- -n * theta^2 / 2 + .evx_log_i2_negative(-z) - chi_constant
  } else {
    q <- .evx_gamma_quad(n, order)
    a <- -z
    rate <- pmax(a, sqrt(n))
    # Q ~ Gamma(n,1); the positive radial integral is evaluated on its
    # natural scale, including when z is very negative.
    val <- outer(1 - a / rate, q$x) -
      outer(1 / (2 * rate^2), q$x^2)
    val <- sweep(val, 2L, q$log_w, "+")
    maximum <- apply(val, 1L, max)
    lmean <- maximum + log(rowSums(exp(val - maximum)))
    ans <- -n * theta^2 / 2 + lgamma(n) - n * log(rate) +
      lmean - chi_constant
  }
  ans[theta == 0] <- 0
  if (anyNA(ans) || any(ans == Inf)) {
    stop("The invariant likelihood ratio exceeded numerical range.", call. = FALSE)
  }
  ans
}

.evx_half_grid <- function(upper, u, n, d, location, scale) {
  like_scale <- if (u > 0) 1 / sqrt(d) else 1 / sqrt(n)
  like_center <- if (u > 0) u * sqrt((n - 1) / (n * d)) else 0
  centers <- unique(c(0, max(location, 0), like_center))
  scales <- unique(c(scale, like_scale))
  offsets <- c(-2^seq(6, -4), 0, 2^seq(-4, 6))
  points <- unlist(lapply(centers, function(center) {
    unlist(lapply(scales, function(s) center + s * offsets), use.names = FALSE)
  }), use.names = FALSE)
  # Bridge wide likelihood/prior scale separations geometrically.
  low <- min(scales) / 16
  if (low > 0 && low < upper) {
    count <- max(2L, ceiling(2 * log2(upper / low)) + 1L)
    if (!is.finite(count) || count > 4096L) {
      stop("The prior/likelihood scale separation is too large.", call. = FALSE)
    }
    points <- c(points, exp(seq(log(low), log(upper), length.out = count)))
  }
  sort(unique(c(0, upper, points[is.finite(points) & points > 0 & points < upper])))
}

.evx_half_once <- function(upper, u, n, d, family, location, scale,
                           log_mass, order) {
  breaks <- .evx_half_grid(upper, u, n, d, location, scale)
  half <- diff(breaks) / 2
  mid <- (breaks[-1L] + breaks[-length(breaks)]) / 2
  good <- half > 0
  half <- half[good]
  mid <- mid[good]
  q <- .evx_legendre(order)
  theta <- as.vector(outer(half, q$x) + mid)
  lw <- as.vector(outer(log(half), q$log_w, "+"))
  lp <- if (family == "gaussian") {
    stats::dnorm(theta, mean = location, sd = scale, log = TRUE)
  } else {
    stats::dcauchy(theta, location = location, scale = scale, log = TRUE)
  }
  .evx_logsum(.evx_log_lr(theta, u, n, d, order) + lp - log_mass + lw)
}

.evx_half_mixture <- function(u, n, d, family, location, scale,
                              log_mass, rel.tol = 1e-8) {
  if (!is.finite(log_mass)) {
    stop("The conditioning probability is outside numerical range.", call. = FALSE)
  }
  # This is a normalized half-prior. A genuine likelihood envelope bounds
  # the omitted mixture independently of the prior density's height.
  log_c <- if (u <= 0) 0 else (n / 2) * log(2 * n / d)
  kappa <- if (u <= 0) n / 2 else n * d / (2 * (2 * n - d))
  upper <- max(1, 4 / sqrt(if (u > 0) d else n))
  for (attempt in seq_len(60L)) {
    coarse <- .evx_half_once(upper, u, n, d, family, location,
                             scale, log_mass, 32L)
    fine <- .evx_half_once(upper, u, n, d, family, location,
                           scale, log_mass, 64L)
    discrepancy <- abs(expm1(coarse - fine))
    if (!is.finite(fine)) {
      stop("The mixture integral is outside numerical range.", call. = FALSE)
    }
    if (!is.finite(discrepancy) || discrepancy > rel.tol / 4) {
      refined <- .evx_half_once(upper, u, n, d, family, location,
                                scale, log_mass, 128L)
      discrepancy <- abs(expm1(fine - refined))
      fine <- refined
      if (!is.finite(discrepancy) || discrepancy > rel.tol / 4) {
        stop("Mixture quadrature did not converge under order refinement.",
             call. = FALSE)
      }
    }
    log_tail <- log_c - kappa * upper^2
    if (log_tail - fine <= log(rel.tol / 4)) return(fine)
    required <- sqrt(max(0, (log_c - fine - log(rel.tol / 8)) / kappa))
    upper <- max(1.5 * upper, 1.05 * required)
    if (!is.finite(upper)) break
  }
  stop("The mixture integral failed its omitted-tail bound.", call. = FALSE)
}

.evx_components <- function(u, n, d, prior, rel.tol = 1e-8) {
  if (!is.finite(u) || !is.finite(d) || d <= 0 || d > n * (1 + 1e-12)) {
    stop("Supply finite invariant information with positive residual energy.",
         call. = FALSE)
  }
  if (prior$family == "lebesgue") {
    stop("Use the separate Lebesgue formulas.", call. = FALSE)
  }
  ans <- c(positive = -Inf, negative = -Inf)
  if (prior$alternative != "less") {
    val <- .evx_half_mixture(u, n, d, prior$family, prior$location,
                             prior$scale, prior$log_positive_mass, rel.tol)
    ans["positive"] <- val + if (prior$alternative == "two.sided") {
      prior$log_positive_mass
    } else 0
  }
  if (prior$alternative != "greater") {
    val <- .evx_half_mixture(-u, n, d, prior$family, -prior$location,
                             prior$scale, prior$log_negative_mass, rel.tol)
    ans["negative"] <- val + if (prior$alternative == "two.sided") {
      prior$log_negative_mass
    } else 0
  }
  ans
}

.evx_log_mixture_stats <- function(st, prior, rel.tol = 1e-8) {
  # d is computed from residuals by .ev_stats, not by subtracting n-U^2.
  h <- max(abs(st$t), sqrt(st$n - 1))
  u <- sqrt(st$n) * (st$t / h) /
    sqrt((st$t / h)^2 + (sqrt(st$n - 1) / h)^2)
  .evx_logsum(.evx_components(u, st$n, st$d, prior, rel.tol))
}

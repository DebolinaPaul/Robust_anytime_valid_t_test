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
  lp <- .evx_prior_logpdf(theta, family, location, scale)
  if (family == "horseshoe") {
    # Apply theta = location +/- r^4 on every interval. The grid includes
    # the location whenever it lies inside the integration range, so no
    # interval crosses the density spike. Retaining the offset separately
    # avoids cancellation when evaluating the density near its location.
    left <- breaks[-length(breaks)][good]
    right <- breaks[-1L][good]
    for (j in seq_along(mid)) {
      direction <- if (left[j] >= location) 1 else -1
      distances <- sort(abs(c(left[j], right[j]) - location))
      r_bounds <- sqrt(sqrt(distances))
      r_half <- (r_bounds[2L] - r_bounds[1L]) / 2
      # Retain the original smooth-interval rule if transformed endpoints
      # cannot be distinguished at floating-point precision.
      if (r_half <= 0) next
      indices <- j + length(mid) * (seq_along(q$x) - 1L)
      r <- r_bounds[1L] + r_half * (1 + q$x)
      offset <- direction * r^4
      theta[indices] <- location + offset
      lw[indices] <- q$log_w + log(r_half) + log(4) + 3 * log(r)
      lp[indices] <- .evx_prior_logpdf(offset, family, 0, scale)
    }
  }
  .evx_logsum(.evx_log_lr(theta, u, n, d, order) + lp - log_mass + lw)
}

.evx_half_mixture <- function(u, n, d, family, location, scale,
                              log_mass, rel.tol = 1e-8) {
  if (!is.finite(log_mass)) {
    stop("The conditioning probability is outside numerical range.", call. = FALSE)
  }
  # This is a normalized half-prior. A likelihood envelope bounds the
  # omitted mixture independently of the prior density's height.
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
    for (order in c(128L, 256L, 512L)) {
      if (is.finite(discrepancy) && discrepancy <= rel.tol / 4) break
      refined <- .evx_half_once(upper, u, n, d, family, location,
                                scale, log_mass, order)
      if (!is.finite(refined)) {
        stop("The mixture integral is outside numerical range.", call. = FALSE)
      }
      discrepancy <- abs(expm1(fine - refined))
      fine <- refined
    }
    if (!is.finite(discrepancy) || discrepancy > rel.tol / 4) {
      stop("Mixture quadrature did not converge under order refinement.",
           call. = FALSE)
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

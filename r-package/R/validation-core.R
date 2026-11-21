# Initial validation core: one-sample, two-sided, centered priors only.
# This file deliberately does not implement source fitting or shifted priors.

.ev_positive_scale <- function(scale, prior) {
  if (is.null(scale)) scale <- if (prior == "rouder") 1 / sqrt(2) else 1
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) ||
      scale <= 0 || !is.finite(scale^2) || scale^2 == 0) {
    stop("`scale` must be positive, with a representable squared value.",
         call. = FALSE)
  }
  scale
}

.ev_logadd <- function(a, b) {
  m <- pmax(a, b)
  ans <- m + log1p(exp(-abs(a - b)))
  ans[is.infinite(m) & m < 0] <- -Inf
  ans
}

.ev_stats <- function(x, mu = 0) {
  if (!is.numeric(x) || !is.null(dim(x)) || length(x) < 2L ||
      any(!is.finite(x)) || !is.numeric(mu) || length(mu) != 1L ||
      !is.finite(mu)) {
    stop("Supply a finite numeric vector of length >= 2 and a finite scalar `mu`.",
         call. = FALSE)
  }
  y <- x - mu
  if (any(!is.finite(y))) stop("Centering overflowed.", call. = FALSE)
  a <- max(abs(y))
  if (a == 0) stop("Positive sample variance is required.", call. = FALSE)
  z <- y / a
  n <- length(z)
  q <- sum(z^2)
  residual <- sum((z - mean(z))^2)
  if (residual <= 0) stop("Positive sample variance is required.", call. = FALSE)
  # This avoids cancellation in n - (sum(z)/sqrt(q))^2.
  d <- n * residual / q
  if (!is.finite(d) || d <= 0 || d > n * (1 + 1e-12)) {
    stop("The standardized residual energy is outside numerical range.",
         call. = FALSE)
  }
  list(n = n, d = min(d, n),
       t = sqrt(n) * mean(z) / sqrt(residual / (n - 1)))
}

.ev_log_gaussian <- function(n, d, scale) {
  ls2 <- 2 * log(scale)
  ((n - 1) / 2) * .ev_logadd(0, log(n) + ls2) -
    (n / 2) * .ev_logadd(0, log(d) + ls2)
}

.ev_log_lebesgue <- function(n, d) {
  0.5 * log(2 * pi / n) + (n / 2) * (log(n) - log(d))
}

.ev_log_rouder <- function(n, d, scale, rel.tol = 1e-10) {
  # theta | W ~ N(0, scale^2/W), W ~ chi-square_1.
  # After cancellation of the sqrt(W) terms, integrate over W >= 0.
  la <- log(n) + 2 * log(scale)
  lb <- log(d) + 2 * log(scale)
  log_h <- function(w) {
    -w / 2 + ((n - 1) / 2) * .ev_logadd(log(w), la) -
      (n / 2) * .ev_logadd(log(w), lb) - 0.5 * log(2 * pi)
  }
  # W = exp(lb) * expm1(u) resolves narrow peaks at W = 0.
  log_f <- function(u) {
    lw <- lb + u + log(-expm1(-u))
    w <- exp(lw)
    -w / 2 + ((n - 1) / 2) * .ev_logadd(lw, la) +
      (1 - n / 2) * (lb + u) - 0.5 * log(2 * pi)
  }
  upper_w <- 64
  for (attempt in seq_len(8L)) {
    upper_u <- .ev_logadd(0, log(upper_w) - lb)
    grid <- seq(0, upper_u, length.out = 65L)
    vals <- log_f(grid)
    j <- which.max(vals)
    lo <- grid[max(1L, j - 1L)]
    hi <- grid[min(length(grid), j + 1L)]
    opt <- stats::optimize(log_f, c(lo, hi), maximum = TRUE)
    shift <- max(vals, opt$objective)
    # Split at the resolved peak; adaptive quadrature sees both sides.
    cuts <- sort(unique(c(0, opt$maximum, upper_u)))
    pieces <- lapply(seq_len(length(cuts) - 1L), function(k) {
      stats::integrate(function(u) exp(log_f(u) - shift),
        cuts[k], cuts[k + 1L], rel.tol = rel.tol,
        abs.tol = 0, subdivisions = 500L, stop.on.error = TRUE)
    })
    value <- sum(vapply(pieces, function(z) z$value, numeric(1)))
    error <- sum(vapply(pieces, function(z) z$abs.error, numeric(1)))
    if (!is.finite(value) || value <= 0 || error > 10 * rel.tol * value) {
      stop("Rouder quadrature failed its numerical error check.", call. = FALSE)
    }
    answer <- shift + log(value)
    # d <= n implies h'(w)/h(w) <= -1/2, hence this tail bound.
    log_tail_bound <- log(2) + log_h(upper_w)
    if (log_tail_bound - answer < log(rel.tol)) return(answer)
    upper_w <- 2 * upper_w
  }
  stop("Rouder integral failed its tail-bound check.", call. = FALSE)
}

#' Log mixture evidence for centered two-sided priors
#'
#' @param x Ordered observations; a finite numeric vector.
#' @param prior One of "gaussian", "rouder", or "lebesgue".
#' @param scale Prior standard deviation (Gaussian) or Cauchy scale (Rouder).
#'   Defaults to 1 and 1/sqrt(2), respectively. Omit for Lebesgue.
#' @param mu Null mean, default zero.
#' @return The natural logarithm of the mixture evidence.
#' @details Only centered, two-sided priors are supported. The effect is
#'   (mean-mu)/sd. Positive sample variance is required. Lebesgue measure is
#'   unnormalized; its evidence must use its own calibration.
#' @seealso [ev_log_mixture()], [ev_ttest()]
#' @export
ev_log_evidence <- function(x, prior = c("gaussian", "rouder", "lebesgue"),
                            scale = NULL, mu = 0) {
  prior <- match.arg(prior)
  st <- .ev_stats(x, mu)
  if (prior == "lebesgue") {
    if (!is.null(scale)) stop("Omit `scale` for Lebesgue measure.", call. = FALSE)
    return(.ev_log_lebesgue(st$n, st$d))
  }
  scale <- .ev_positive_scale(scale, prior)
  if (prior == "gaussian") return(.ev_log_gaussian(st$n, st$d, scale))
  .ev_log_rouder(st$n, st$d, scale)
}

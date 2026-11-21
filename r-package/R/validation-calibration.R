# Calibration from n0 = 2. The Gaussian and Lebesgue expressions are analytic.
# Rouder calibration uses adaptive quadrature, not a certified interval bound.

.ev_g_gaussian <- function(log_c, scale) {
  b <- 2 * scale^2
  if (!is.finite(b)) stop("Gaussian calibration scale is too large.", call. = FALSE)
  log_a <- 0.5 * log1p(b)
  ans <- rep(1, length(log_c))
  high <- log_c >= log_a
  ans[high] <- exp(-log_c[high])
  mid <- log_c > -log_a & log_c < log_a
  if (any(mid)) {
    lc <- log_c[mid]
    q <- pmin(1, pmax(0, expm1(log_a - lc) / b))
    cosine <- sqrt(q)
    sine <- sqrt(pmax(0, 1 - cosine^2))
    v <- atan2(sine, cosine)
    primitive <- atan2(sine, exp(log_a) * cosine)
    ans[mid] <- (2 / pi) * (exp(-lc) * primitive + asin(cosine))
  }
  ans
}

.ev_g_lebesgue_log <- function(log_c) {
  # Form sqrt(q) directly: q itself can underflow while g(c) is representable.
  z <- exp(pmin(0, (0.5 * log(pi) - log_c) / 2))
  (2 / pi) * (asin(z) + z * sqrt(pmax(0, 1 - z^2)))
}

# Stable log(atan(exp(z))), including very small angles.
.ev_log_atan_exp <- function(z) {
  ans <- z
  mid <- is.finite(z) & z >= -20 & z < 0
  ans[mid] <- log(atan(exp(z[mid])))
  high <- z >= 0
  ans[high] <- log(atan2(1, exp(-z[high])))
  ans
}

.ev_g_rouder <- function(log_c, scale) {
  minimum <- .ev_log_rouder(2, 2, scale)
  vapply(log_c, function(lc) {
    if (lc <= minimum) return(1)
    if (lc == Inf) return(0)
    # Solve M_2 = c in log residual energy, not near an angular endpoint.
    objective <- function(ld) .ev_log_rouder(2, exp(ld), scale) - lc
    lower <- -2
    while (objective(lower) < 0 && lower > -700) {
      lower <- max(-700, 2 * lower)
    }
    if (objective(lower) < 0) {
      # The crossing cannot be represented at this precision. The ordinary
      # Ville bound remains a valid (conservative) calibration.
      warning("Rouder crossing unresolved; using the conservative Ville bound.",
              call. = FALSE)
      return(exp(min(0, -lc)))
    }
    ld <- stats::uniroot(objective, c(lower, log(2)), tol = 1e-11)$root
    q <- min(1, exp(ld - log(2)))
    delta <- asin(sqrt(q))
    if (q >= 1) return(1)
    log_sin <- 0.5 * (ld - log(2))
    log_cos <- 0.5 * log1p(-q)
    # Integrate the truncated first moment analytically in the null angle,
    # then numerically over the half-normal scale variable Z.
    log_phi <- function(z) {
      lz <- log(z)
      log_ratio <- lz + log_cos - log_sin -
        0.5 * .ev_logadd(2 * lz, log(2) + 2 * log(scale))
      .ev_log_atan_exp(log_ratio)
    }
    shift <- log_phi(1)
    out <- stats::integrate(function(z) {
      exp(-z^2 / 2 + log_phi(z) - shift)
    }, 0, Inf, rel.tol = 1e-10, abs.tol = 0,
       subdivisions = 500L, stop.on.error = TRUE)
    if (!is.finite(out$value) || out$value <= 0 ||
        out$abs.error > 1e-9 * out$value) {
      stop("Rouder calibration failed its relative-error check.", call. = FALSE)
    }
    log_partial_mean <- 0.5 * log(2 / pi) + shift + log(out$value)
    exp(.ev_logadd(log(2 / pi) + log(delta),
                  log(2 / pi) - lc + log_partial_mean))
  }, numeric(1))
}

#' Calibrate two-sided evidence starting at n = 2
#'
#' @param log_evidence Natural logarithms of evidence, possibly infinite.
#' @param prior One of "gaussian", "rouder", or "lebesgue".
#' @param scale Prior scale, as in ev_log_evidence().
#' @param method Either "eVille" or "Ville".
#' @return Calibrated values between 0 and 1, up to numerical roundoff.
#' @details Proper priors use the null law corresponding to that same prior.
#'   Rouder calibration is numerical. Ordinary Ville is unavailable for
#'   Lebesgue measure. Pass the cumulative maximum of log evidence to obtain
#'   running p-values; pass the current log evidence for pointwise p-values.
#' @export
ev_calibrate <- function(log_evidence,
                         prior = c("gaussian", "rouder", "lebesgue"),
                         scale = NULL, method = c("eVille", "Ville")) {
  prior <- match.arg(prior)
  method <- match.arg(method)
  if (!is.numeric(log_evidence) || anyNA(log_evidence)) {
    stop("`log_evidence` must be numeric without NA or NaN.", call. = FALSE)
  }
  if (prior == "lebesgue") {
    if (!is.null(scale)) stop("Omit `scale` for Lebesgue measure.", call. = FALSE)
    if (method == "Ville") {
      stop("Ordinary Ville calibration is not defined for Lebesgue measure.",
           call. = FALSE)
    }
    return(.ev_g_lebesgue_log(log_evidence))
  }
  scale <- .ev_positive_scale(scale, prior)
  if (method == "Ville") return(exp(pmin(0, -log_evidence)))
  if (prior == "gaussian") return(.ev_g_gaussian(log_evidence, scale))
  .ev_g_rouder(log_evidence, scale)
}

#' Evidence and running p-values for a two-sided one-sample test
#'
#' @param x Ordered observations. Every prefix from n = 2 must have positive
#'   sample variance.
#' @param prior Prior family, as in ev_log_evidence().
#' @param scale Prior scale, as in ev_log_evidence().
#' @param mu Null mean.
#' @return A data frame with n, log_evidence, p_eville, and p_ville.
#' @details Monitoring begins at 2; order is preserved. For Lebesgue, p_ville
#'   is NA. This helper supports centered two-sided priors. Use
#'   [ev_ttest()] for fitted objects, shifted priors, or one-sided alternatives.
#' @seealso [ev_ttest()], [ev_log_mixture()]
#' @export
ev_path <- function(x, prior = c("gaussian", "rouder", "lebesgue"),
                     scale = NULL, mu = 0) {
  prior <- match.arg(prior)
  .ev_stats(x, mu)
  n <- seq.int(2L, length(x))
  lm <- vapply(n, function(k) {
    ev_log_evidence(x[seq_len(k)], prior, scale, mu)
  }, numeric(1))
  running <- cummax(lm)
  pe <- ev_calibrate(running, prior, scale)
  pv <- if (prior == "lebesgue") rep(NA_real_, length(n)) else
    ev_calibrate(running, prior, scale, method = "Ville")
  data.frame(n = n, log_evidence = lm, p_eville = pe, p_ville = pv)
}

.evx_log_half_flat <- function(st, direction = 1) {
  n <- st$n
  t <- direction * st$t
  if (n == 2) {
    return(0.5 * log(pi) - log(4) + .ev_logadd(0, 2 * asinh(t)))
  }
  if (t >= -10) {
    return(.ev_log_lebesgue(n, st$d) +
      stats::pt(t * sqrt(n / (n - 1)), df = n, log.p = TRUE))
  }
  # Cancel the opposing large powers before evaluating the negative tail.
  z <- st$d / n
  a <- n / 2
  h <- stats::integrate(function(q) {
    exp(-q) / sqrt((1 - z) + z * (-expm1(-q / a)))
  }, 0, Inf, rel.tol = 1e-10, abs.tol = 0,
     subdivisions = 500L, stop.on.error = TRUE)
  if (h$abs.error > 1e-9 * h$value) stop("Half-Lebesgue quadrature failed.")
  0.5 * log(2 * pi / n) - log(n) - lbeta(n / 2, 0.5) + log(h$value)
}

.evx_log_mixture_dispatch <- function(st, prior) {
  if (!prior$proper) {
    if (prior$alternative == "two.sided") return(.ev_log_lebesgue(st$n, st$d))
    return(.evx_log_half_flat(st, if (prior$alternative == "greater") 1 else -1))
  }
  if (prior$alternative == "two.sided" && prior$location == 0) {
    if (prior$family == "gaussian") {
      return(.ev_log_gaussian(st$n, st$d, prior$scale))
    }
    return(.ev_log_rouder(st$n, st$d, prior$scale))
  }
  .evx_log_mixture_stats(st, prior)
}

#' Compute evidence under a shifted or one-sided prior
#'
#' @param x Finite ordered observations with positive sample variance.
#' @param prior A fixed prior constructed by ev_prior().
#' @param mu Null mean, default zero.
#' @return Natural logarithm of the mixture evidence.
#' @export
ev_log_mixture <- function(x, prior = ev_prior(), mu = 0) {
  .evx_check_prior(prior)
  .evx_log_mixture_dispatch(.ev_stats(x, mu), prior)
}

.evx_first_rejection <- function(n, p, alpha) {
  i <- which(!is.na(p) & p <= alpha)
  if (length(i)) n[i[1L]] else NA_integer_
}

#' Run an anytime-valid one-sample t-test
#'
#' @param x Ordered finite observations. Each prefix from n=2 must have
#'   positive sample variance; observations are never silently removed.
#' @param prior A fixed prior constructed with ev_prior().
#' @param mu Null boundary mean.
#' @param alpha Significance level between zero and one.
#' @param n0 Initial monitoring time; this implementation requires 2.
#' @param calibration Optional reusable object from ev_calibration(), made
#'   for exactly the same prior, support, and monitoring start.
#' @param tolerance Absolute tolerance for the null-angle calibration bracket
#'   when constructing a new calibrator. Defaults to 1e-4.
#' @return An ev_ttest object with path, stopping_times, prior, calibration,
#'   alpha, mu, and call. Ville is NA for Lebesgue measure.
#' @details The working model is independent normal observations with unknown
#'   common variance. Two-sided tests use the point null mu=mu0. Greater tests
#'   use null mu<=mu0, and less tests use null mu>=mu0; zero standardized effect
#'   is the least-favourable boundary for either one-sided null. General
#'   calibration returns a conservative upper approximation, together with
#'   its null-angle integration gap. Floating-point quadrature error in the
#'   evidence is separate from that gap. Source-fitted priors require source
#'   data independent of the target, and are held fixed while monitoring.
#'   A warning identifies prefixes where the integration bracket straddles
#'   alpha; reduce tolerance to resolve a decision at such a prefix.
#' @seealso [ev_prior()], [ev_calibration()], [ev_log_mixture()]
#' @examples
#' x <- c(-0.3, -1.1, -0.7, -0.5, -1.4, -0.8, -0.2, -1.0)
#' prior <- ev_prior("gaussian", location = -4, scale = 1,
#'                   alternative = "less")
#' fit <- ev_ttest(x, prior = prior, alpha = 0.05)
#' summary(fit)
#' plot(fit)
#' head(as.data.frame(fit))
#' @export
ev_ttest <- function(x, prior = ev_prior(), mu = 0, alpha = 0.05, n0 = 2,
                     calibration = NULL, tolerance = 1e-4) {
  .evx_check_prior(prior)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1 || length(n0) != 1L || is.na(n0) || n0 != 2) {
    stop("Require 0 < alpha < 1 and n0 = 2.", call. = FALSE)
  }
  .ev_stats(x, mu)
  if (is.null(calibration)) calibration <- ev_calibration(prior, tolerance)
  if (!inherits(calibration, "ev_calibration") ||
      !identical(calibration$prior, prior)) {
    stop("The calibration must match this exact prior and alternative.",
         call. = FALSE)
  }
  n <- seq.int(2L, length(x))
  lm <- vapply(n, function(k) {
    .evx_log_mixture_dispatch(.ev_stats(x[seq_len(k)], mu), prior)
  }, numeric(1))
  running <- cummax(lm)
  cal <- ev_pvalues(running, calibration)
  pv <- if (prior$proper) exp(pmin(0, -running)) else rep(NA_real_, length(n))
  path <- data.frame(n = n, log_evidence = lm,
    p_eville = cal$p_eville, p_ville = pv,
    calibration_lower = cal$lower, calibration_gap = cal$gap)
  ambiguous <- n[cal$lower <= alpha & cal$p_eville > alpha]
  if (length(ambiguous)) {
    warning(paste0("The calibration bracket straddles alpha at n = ",
      paste(ambiguous[seq_len(min(8L, length(ambiguous)))], collapse = ", "),
      if (length(ambiguous) > 8L) ", ..." else "",
      ". The test uses the conservative upper value; use a tighter calibration to resolve it."),
      call. = FALSE)
  }
  stopping <- c(eVille = .evx_first_rejection(n, path$p_eville, alpha),
                 Ville = .evx_first_rejection(n, path$p_ville, alpha))
  structure(list(call = match.call(), prior = prior, calibration = calibration,
    mu = mu, alpha = alpha, path = path, stopping_times = stopping,
    calibration_ambiguous_n = ambiguous),
    class = "ev_ttest")
}

#' @export
print.ev_ttest <- function(x, ...) {
  cat("Anytime-valid one-sample t-test\n")
  cat("Alternative:", x$prior$alternative, "| Null boundary:", x$mu, "\n")
  cat("Prior:", x$prior$family)
  if (x$prior$proper) {
    cat(" | location:", x$prior$location, "| scale:", x$prior$scale)
  }
  cat("\nLevel:", x$alpha, "| Monitoring starts at n = 2\n")
  print(x$stopping_times)
  if (x$prior$proper) cat("NA means no rejection in the observed sequence.\n")
  else cat("Ville is unavailable for Lebesgue; eVille NA means no rejection.\n")
  invisible(x)
}

#' @export
as.data.frame.ev_ttest <- function(x, ...) x$path

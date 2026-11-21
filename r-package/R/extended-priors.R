#' Specify a fixed effect-size prior
#'
#' @param family One of "gaussian", "cauchy", "rouder", or "lebesgue".
#' @param location Location of the unconditioned distribution on the
#'   standardized effect (mu-mu0)/sigma. Rouder requires location zero.
#' @param scale Gaussian standard deviation or Cauchy scale; defaults to
#'   1 for Gaussian and 1/sqrt(2) for Cauchy/Rouder. Omit for Lebesgue.
#' @param alternative One of "two.sided", "greater", or "less".
#' @return An object of class ev_prior describing the fixed mixing measure.
#' @details One-sided proper priors are the specified distribution
#'   conditioned on the selected half-line. Location and scale describe the
#'   distribution before conditioning, not its conditional mean and SD.
#'   Fix the prior before examining target observations. Lebesgue measure is
#'   restricted to the selected half-line without normalization.
#' @export
ev_prior <- function(family = c("gaussian", "cauchy", "rouder", "lebesgue"),
                     location = 0, scale = NULL,
                     alternative = c("two.sided", "greater", "less")) {
  family <- match.arg(family)
  alternative <- match.arg(alternative)
  if (!is.numeric(location) || length(location) != 1L || !is.finite(location)) {
    stop("`location` must be one finite number.", call. = FALSE)
  }
  if (family == "lebesgue") {
    if (location != 0 || !is.null(scale)) {
      stop("Omit location and scale for Lebesgue measure.", call. = FALSE)
    }
    return(structure(list(family = family, location = 0, scale = NULL,
      alternative = alternative, proper = FALSE,
      log_positive_mass = NA_real_, log_negative_mass = NA_real_),
      class = "ev_prior"))
  }
  if (family == "rouder" && location != 0) {
    stop("For a shifted Cauchy use family = 'cauchy'.", call. = FALSE)
  }
  if (family == "rouder") family <- "cauchy"
  scale <- .ev_positive_scale(scale, if (family == "cauchy") "rouder" else family)
  if (location != 0 && scale <= 16 * .Machine$double.eps * abs(location)) {
    stop("The prior scale is below numerical resolution at this location.",
         call. = FALSE)
  }
  if (!is.finite(location / scale)) {
    stop("The location/scale ratio exceeds numerical range.", call. = FALSE)
  }
  cdf <- if (family == "gaussian") stats::pnorm else stats::pcauchy
  positive <- cdf(-location / scale, lower.tail = FALSE, log.p = TRUE)
  negative <- cdf(-location / scale, lower.tail = TRUE, log.p = TRUE)
  if (any(!is.finite(c(positive, negative)))) {
    stop("A prior half-line mass is outside numerical range.", call. = FALSE)
  }
  structure(list(family = family, location = location, scale = scale,
    alternative = alternative, proper = TRUE,
    log_positive_mass = positive, log_negative_mass = negative),
    class = "ev_prior")
}

.evx_check_prior <- function(prior) {
  if (!inherits(prior, "ev_prior")) {
    stop("Construct the prior with ev_prior().", call. = FALSE)
  }
  rebuilt <- ev_prior(prior$family, prior$location, prior$scale,
                       prior$alternative)
  if (!identical(prior, rebuilt)) {
    stop("The prior object was modified; construct a fresh ev_prior().",
         call. = FALSE)
  }
  invisible(prior)
}

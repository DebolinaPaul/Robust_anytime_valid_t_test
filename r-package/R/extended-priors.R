#' Specify a fixed effect-size prior
#'
#' @param family One of "gaussian", "cauchy", "student_t5", "laplace",
#'   "horseshoe", "rouder", or "lebesgue". The aliases "student_t",
#'   "student-t5", and "t5" select the Student family with five degrees of
#'   freedom. Names are case insensitive.
#' @param location Location of the unconditioned distribution on the
#'   standardized effect (mu-mu0)/sigma. Rouder requires location zero.
#' @param scale Scale of the unconditioned family: Gaussian standard
#'   deviation, Cauchy scale, Student-t scale, Laplace scale, or horseshoe
#'   global scale. Defaults to 1/sqrt(2) for Cauchy/Rouder and 1 otherwise,
#'   matching the simulation notebook. Omit for Lebesgue.
#' @param alternative One of "two.sided", "greater", or "less".
#' @return An object of class ev_prior describing the fixed mixing measure.
#' @details One-sided proper priors are the specified distribution
#'   conditioned on the selected half-line. Location and scale describe the
#'   distribution before conditioning, not its conditional mean and SD.
#'   Fix the prior before examining target observations. Lebesgue measure is
#'   restricted to the selected half-line without normalization.
#'   Student-t5 has exactly five degrees of freedom; its scale is not its
#'   standard deviation. Laplace has density exp(-abs(theta-location)/scale)
#'   divided by twice scale. The genuine horseshoe is the normal-scale
#'   mixture theta given lambda distributed as N(location,scale^2*lambda^2),
#'   with lambda half-Cauchy(0,1); its density is infinite at its location.
#' @export
ev_prior <- function(family = c("gaussian", "cauchy", "student_t5", "laplace",
                                "horseshoe", "rouder", "lebesgue"),
                     location = 0, scale = NULL,
                     alternative = c("two.sided", "greater", "less")) {
  if (missing(family)) family <- "gaussian"
  if (!is.character(family) || length(family) != 1L || is.na(family)) {
    stop("`family` must be one prior-family name.", call. = FALSE)
  }
  family <- tolower(family)
  if (family %in% c("student_t", "student-t5", "t5")) family <- "student_t5"
  family <- match.arg(family, c("gaussian", "cauchy", "student_t5", "laplace",
                                "horseshoe", "rouder", "lebesgue"))
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
  masses <- .evx_prior_log_masses(family, location, scale)
  positive <- unname(masses["positive"])
  negative <- unname(masses["negative"])
  if (any(!is.finite(c(positive, negative)))) {
    stop("A prior half-line mass is outside numerical range.", call. = FALSE)
  }
  structure(list(family = family, location = location, scale = scale,
    alternative = alternative, proper = TRUE,
    log_positive_mass = positive, log_negative_mass = negative),
    class = "ev_prior")
}

#' Evaluate the specified effect-size prior density
#'
#' @param theta Finite numeric effect sizes.
#' @param prior A fixed prior constructed with [ev_prior()].
#' @param log Return the natural logarithm of the density or intensity.
#' @return A numeric vector. Proper priors return their normalized density,
#'   including the normalization for a one-sided restriction. Lebesgue
#'   measure returns intensity one on its support and zero off its support;
#'   it is not a probability density. The genuine horseshoe density is
#'   infinite at its location when that point is in the support.
#' @seealso [ev_prior()], [ev_log_mixture()]
#' @export
ev_prior_density <- function(theta, prior = ev_prior(), log = FALSE) {
  .evx_check_prior(prior)
  if (!is.numeric(theta) || any(!is.finite(theta)) ||
      !is.logical(log) || length(log) != 1L || is.na(log)) {
    stop("Supply finite numeric effects and one logical `log` value.", call. = FALSE)
  }
  inside <- rep(TRUE, length(theta))
  mass <- 0
  if (prior$alternative == "greater") {
    inside <- theta > 0
    mass <- prior$log_positive_mass
  } else if (prior$alternative == "less") {
    inside <- theta < 0
    mass <- prior$log_negative_mass
  }
  answer <- rep(-Inf, length(theta))
  if (!prior$proper) answer[inside] <- 0
  else if (any(inside)) {
    answer[inside] <- .evx_prior_logpdf(theta[inside], prior$family,
      prior$location, prior$scale) - mass
  }
  if (log) answer else exp(answer)
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

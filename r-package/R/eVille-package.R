#' eVille: Anytime-Valid One-Sample t-Tests
#'
#' @description
#' One-sample t-tests for a Gaussian mean with unknown variance, using fixed
#' effect-size priors and Ville or extended Ville calibration. The main
#' interface is [ev_ttest()], with [ev_prior()] specifying the prior and
#' alternative. Proper priors use Monte Carlo calibration by default;
#' this is a finite-sample approximation to the theoretically anytime-valid
#' population calibration. Use summary() and plot() to inspect the results.
#'
#' @section Model and hypotheses:
#' Observations follow the working model
#' \eqn{Y_i \stackrel{\mathrm{iid}}{\sim} N(\mu,\sigma^2)}, with
#' \eqn{\sigma>0} unknown. The argument `mu` specifies the null boundary
#' \eqn{\mu_0}. The prior is on the standardized effect
#' \eqn{\theta=(\mu-\mu_0)/\sigma}.
#'
#' The alternative in [ev_prior()] determines the hypotheses:
#' \describe{
#'   \item{`two.sided`}{\eqn{H_0:\mu=\mu_0} against \eqn{H_1:\mu\ne\mu_0}.}
#'   \item{`less`}{\eqn{H_0:\mu\ge\mu_0} against \eqn{H_1:\mu<\mu_0}.}
#'   \item{`greater`}{\eqn{H_0:\mu\le\mu_0} against \eqn{H_1:\mu>\mu_0}.}
#' }
#'
#' @section Prior families:
#' Gaussian, Cauchy, Student-t5, Laplace and horseshoe priors use the same
#' location-scale parameterizations as the accompanying Python simulation
#' notebook. The horseshoe is the genuine normal--half-Cauchy scale mixture.
#' The centered Rouder default is a Cauchy with scale \eqn{1/\sqrt{2}}.
#' One-sided proper priors are conditioned on the selected half-line.
#' Location and scale describe the full-line distribution before conditioning.
#' Use [ev_prior_density()] to evaluate the specified density.
#'
#' @section Sequential use:
#' Supply observations in their sampling order. Monitoring starts at `n0 = 2`.
#' Under the stated model, the population calibration controls rejection over
#' all monitoring times. Finite Monte Carlo calibration approximates that
#' construction; it does not itself establish an exact anytime-valid guarantee.
#' The prior and direction must be fixed before examining the target data.
#' Parameters fitted using independent source data can be supplied to
#' [ev_prior()]. Reordering an observed dataset to obtain a smaller value
#' is not covered by the fixed-order guarantee.
#'
#' @section Calibration and numerical precision:
#' The default proper-prior calibration uses the same 100,000 null t
#' statistics as the Python notebook's paper configuration. [ev_null_sample()]
#' returns those bundled draws. A calibrator can be reused for the same
#' prior, and different priors can share the same null sample.
#'
#' The empirical clipped mean is evaluated with numerical lower and upper
#' bounds. The upper value is capped by the ordinary Ville bound. The
#' empirical numerical gap and the Monte Carlo standard error describe
#' different uncertainties: neither is a simultaneous confidence bound
#' correcting the Monte Carlo approximation. Evidence quadrature error is
#' separate from both diagnostics. Smaller numerical tolerance does not
#' reduce Monte Carlo sampling error.
#'
#' Select `calibration_method = "deterministic"` in [ev_ttest()] to use
#' analytic formulas or numerical bounds for the population null expectation.
#' Numerical bounds assume accurate evidence evaluations and are not
#' interval-arithmetic certificates. Lebesgue benchmarks always use their
#' own analytic calibration; ordinary Ville is unavailable for them.
#'
#' @section Output:
#' An [ev_ttest()] result contains running values, mixture log-evidence,
#' first rejection times, the calibration method and diagnostics.
#' `summary(fit)` gives a compact table, `plot(fit)` draws running values,
#' and `as.data.frame(fit)` returns the full trajectory.
#' With a proper prior, eVille and Ville use the same evidence. For Monte
#' Carlo calibration, the reported eVille values are approximations.
#'
#' @section Current scope:
#' The implementation supports the five proper families and their one-sided
#' restrictions, plus full and half-line Lebesgue benchmarks. Every monitored
#' prefix must have positive sample variance. Missing and nonfinite
#' observations are rejected.
#' Two-sample and regression tests, other monitoring start times, automatic
#' prior fitting and incremental streaming updates are not implemented.
#' Robustness to prior-data conflict does not remove the Gaussian or
#' independence assumptions.
#'
#' @seealso [ev_ttest()], [ev_prior()], [ev_prior_density()],
#'   [ev_calibration()], [ev_pvalues()], [ev_null_sample()], [ev_log_mixture()]
#'
#' @examples
#' x <- c(-0.3, -1.1, -0.7, -0.5, -1.4, -0.8, -0.2, -1.0)
#' prior <- ev_prior("gaussian", location = -4, scale = 1,
#'                   alternative = "less")
#' fit <- ev_ttest(x, prior = prior, mu = 0, alpha = 0.05)
#' summary(fit)
#' plot(fit)
#'
#' @keywords internal
"_PACKAGE"

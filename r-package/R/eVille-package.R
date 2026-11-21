#' eVille: Anytime-Valid One-Sample t-Tests
#'
#' @description
#' Anytime-valid one-sample t-tests for a Gaussian mean with unknown variance,
#' using fixed effect-size priors and extended Ville calibration. The main
#' interface is [ev_ttest()], with [ev_prior()] specifying the prior and
#' alternative. Use summary() and plot() to inspect the results.
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
#' @section Sequential use:
#' Supply observations in their sampling order. Monitoring starts at `n0 = 2`.
#' The calibration controls rejection over monitoring times under the stated
#' model; it is not repeated use of ordinary fixed-sample t-test p-values.
#' The prior and direction must be fixed before examining the target data.
#' Parameters fitted using independent source data can be supplied to
#' [ev_prior()]. Reordering an observed dataset to obtain a smaller p-value
#' is not covered by the fixed-order guarantee.
#'
#' @section Output:
#' An [ev_ttest()] result contains running p-values, mixture log-evidence,
#' first rejection times, and calibration diagnostics. `summary(fit)` gives
#' a compact table, `plot(fit)` draws the running p-values, and
#' `as.data.frame(fit)` returns the full trajectory.
#' With a proper prior, eVille and Ville use the same evidence. Lebesgue
#' measure has its own extended Ville calibration; ordinary Ville is
#' unavailable for this benchmark.
#'
#' @section Numerical precision:
#' Analytic calibration is used where implemented. General calibration uses
#' an upper approximation and reports a null-integration gap. This gap does
#' not include every floating-point error in computing the evidence. A warning
#' identifies times at which the calibration bracket straddles the requested
#' significance level; refine the tolerance to resolve those decisions.
#'
#' @section Current scope:
#' The implementation supports Gaussian and Cauchy priors (including the
#' centered Rouder default), their one-sided restrictions, and full or half-line
#' Lebesgue benchmarks. Every monitored prefix must have positive sample
#' variance. Missing and nonfinite observations are rejected.
#' Two-sample and regression tests, other monitoring start times, automatic
#' prior fitting, and incremental streaming updates are not implemented.
#' Robustness to prior-data conflict does not remove the Gaussian or
#' independence assumptions.
#'
#' @seealso [ev_ttest()], [ev_prior()], [ev_calibration()], [ev_pvalues()],
#'   [ev_log_mixture()]
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

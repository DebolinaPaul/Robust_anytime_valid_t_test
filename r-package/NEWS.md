# eVille 0.2.0

* Adds Student-t5, Laplace and genuine normal--half-Cauchy horseshoe priors,
  using the same parameterizations as the Python simulation notebook.
* Adds prior-density evaluation and supports one-sided conditioning of all
  five proper prior families.
* Makes Monte Carlo calibration the default for proper priors, using the
  notebook's shared 100,000 null draws and empirical-calibration rule.
* Reports empirical-calibration brackets and Monte Carlo standard errors;
  explicitly distinguishes finite Monte Carlo approximation from the
  theoretical anytime-valid guarantee.
* Retains deterministic calibration as an explicit option and analytic
  calibration for the full and half-line Lebesgue benchmarks.
* Adds all-prior examples and cross-language evidence/calibration references.

# eVille 0.1.0

* Implements anytime-valid one-sample t-tests for a Gaussian mean with
  unknown variance, with monitoring from the second observation.
* Supports fixed Gaussian and Cauchy priors, including shifted priors and
  the centered Rouder default, with one-sided or two-sided alternatives.
* Provides Ville and extended Ville calibration for proper prior mixtures.
* Includes separately calibrated full and half-line Lebesgue benchmarks;
  ordinary Ville is unavailable for these improper measures.
* Reports running p-values, first rejection times, and calibration diagnostics.
* Provides summary(), plot(), and as.data.frame() methods for test results.
* Adds a package overview and a reproducible synthetic example.

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

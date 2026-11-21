# eVille: Anytime-Valid One-Sample t-Tests

**eVille implements one-sample t-tests for a Gaussian mean with unknown
variance**, using fixed effect-size priors and Ville or extended Ville
calibration. Supply observations in their sampling order to obtain running
values, first rejection times, summaries and plots. Monitoring starts at
observation 2.

Version 0.2.0 includes all five prior families used in the Python simulation
notebook: Gaussian, Cauchy, Student-t5, Laplace and horseshoe. It also includes
the centered Rouder default and full or half-line Lebesgue benchmarks.

**Calibration matters:** proper priors use Monte Carlo calibration by default,
matching the simulation notebook. A finite Monte Carlo approximation does
not itself establish the exact anytime-valid guarantee of the population
calibration. The package reports Monte Carlo and numerical diagnostics;
deterministic calibration remains available explicitly.

## Installation

```r
install.packages("remotes")  # Run once if needed.
remotes::install_github(
  "DebolinaPaul/Robust_anytime_valid_t_test",
  subdir = "r-package",
  upgrade = "never"
)
library(eVille)
```

To install a particular archived version, supply `ref` with its Git tag or
commit SHA. To install local source, run from the directory containing
`DESCRIPTION`:

```r
devtools::install(".", upgrade = FALSE)
```

## Example

The observations below are synthetic. The prior is specified in advance,
independently of these observations.

```r
library(eVille)
x <- c(-0.3, -1.1, -0.7, -0.5, -1.4, -0.8, -0.2, -1.0)

prior <- ev_prior("gaussian", location = -4, scale = 1,
                  alternative = "less")
fit <- ev_ttest(x, prior = prior, mu = 0, alpha = 0.05)

summary(fit)
plot(fit)
summary(fit)$results
head(as.data.frame(fit))
```

The null is that the mean is nonnegative, against a negative mean. The
prior is N(-4, 1), conditioned on negative standardized effects. The eVille
values in this example use the default Monte Carlo approximation. The
ordinary Ville comparison uses the same mixture evidence.

## The five priors

The prior is on `theta = (mu - mu0) / sigma`; `mu` in `ev_ttest()` specifies
`mu0`. Let `m = location` and `s = scale`. The parameterizations below are
exactly those used by the Python simulation notebook.

| Family in `ev_prior()` | Distribution before one-sided conditioning | Default scale |
|---|---|---:|
| `"gaussian"` | Normal with mean `m` and standard deviation `s` | 1 |
| `"cauchy"` | Cauchy with location `m` and scale `s` | `1 / sqrt(2)` |
| `"student_t5"` | `m + s * T`, where `T` has Student's t distribution with 5 degrees of freedom | 1 |
| `"laplace"` | Density `exp(-abs(theta - m) / s) / (2 * s)` | 1 |
| `"horseshoe"` | `theta = m + s * lambda * Z`, with independent `Z ~ N(0,1)` and `lambda ~ half-Cauchy(0,1)` | 1 |

The Student-t5 scale is not its standard deviation; its variance is
`5 * s^2 / 3`. The Laplace variance is `2 * s^2`. The horseshoe is the
normal--half-Cauchy mixture above, including its infinite density at `m`.
Its scale is not a standard deviation.

```r
priors <- list(
  Gaussian = ev_prior("gaussian", location = 3),
  Cauchy = ev_prior("cauchy", location = 3),
  Student_t5 = ev_prior("student_t5", location = 3),
  Laplace = ev_prior("laplace", location = 3),
  Horseshoe = ev_prior("horseshoe", location = 3)
)

# Inspect the actual prior density, including one-sided conditioning.
ev_prior_density(c(-1, 0, 1), ev_prior("student_t5"))

# Additional benchmarks.
p_rouder <- ev_prior("rouder")
p_flat <- ev_prior("lebesgue", alternative = "less")
```

`"rouder"` selects a centered Cauchy, with default scale `1 / sqrt(2)`.
`"lebesgue"` specifies an improper mixing measure, not a probability density;
it is calibrated separately using the analytic full or half-line formula.
Ordinary Ville is unavailable for Lebesgue measure.

| `alternative` | Null hypothesis | Alternative |
|---|---|---|
| `"two.sided"` | mu = mu0 | mu != mu0 |
| `"less"` | mu >= mu0 | mu < mu0 |
| `"greater"` | mu <= mu0 | mu > mu0 |

Every proper family supports these alternatives. One-sided priors are the
specified full-line distribution conditioned on the chosen half-line.
Location and scale refer to the distribution before conditioning.
The simulation notebook uses full-line, two-sided priors; one-sided support
is an additional package feature.

## Monte Carlo calibration and agreement with Python

For each proper prior, the ideal extended Ville calibration is

$$
g(c) = \mathbb E_0[\min\{1,M_2/c\}],
\qquad c>0.
$$

The default implementation approximates this expectation with the same
100,000 independent null pairs used in the notebook's paper configuration
(NumPy seed `20260510`). Their derived t statistics are bundled with the
package. Consequently, using the defaults does not depend on R and NumPy
producing the same random sequence from the same integer seed.

For the finite sample of null evidence values, the empirical calibration is

$$
\widehat g_B(c)=\frac1B\sum_{b=1}^{B}\min\{1,M_{2,b}^{(0)}/c\}.
$$

As in the Python notebook, adaptive convex bounds evaluate this finite-sample
mean. The value used is the numerical upper estimate capped by the exact
ordinary Ville bound `min(1, 1/c)`. This preserves `p_eville <= p_ville` for
the same proper mixture. The cap does not remove Monte Carlo error or turn
the empirical procedure into an exactly calibrated population procedure.

```r
null_t <- ev_null_sample()  # The same 100,000 null t statistics as Python.
prior <- ev_prior("student_t5", location = 3)
cal <- ev_calibration(prior, method = "monte_carlo", null_t = null_t)

fit <- ev_ttest(x, prior = prior, calibration = cal)
summary(fit)

# Reuse this calibrator for further target samples with the same prior.
# For another prior, construct its own calibrator using the same null_t.
```

The default absolute calibration tolerance is `1e-5`, matching the notebook's
paper configuration. Running trajectories also use `relative_tolerance = 1e-3`,
matching its trajectory calculation. Both requirements apply to the numerical
empirical-calibration bracket. Standalone `ev_pvalues()` calls default to
absolute tolerance only; supply `relative_tolerance = 1e-3` when reproducing
a trajectory calculation.

Calibration objects cache evidence evaluations and refine their numerical
bounds when needed. Reusing a calibrator avoids repeating those calculations.

A custom vector of independent null t statistics may be supplied as
`null_t`. It must be generated independently of the target data and must
have the null t distribution with one degree of freedom. Such a custom
sample changes the Monte Carlo results; giving R and NumPy the same seed
alone does not reproduce the bundled sample.

Matching requires the same observations, prior parameters, alternative,
null sample and numerical tolerances. Floating-point implementations and
quadrature can still differ within their numerical accuracy. The
cross-language tests compare evidence and empirical calibration against
references exported from the current Python notebook.

## Numerical and Monte Carlo uncertainty

There are three distinct sources of numerical uncertainty:

1. **Evidence quadrature.** Refinement checks and omitted-tail bounds assess
   the mixture calculation. They are not interval-arithmetic certificates.
2. **Empirical-calibration integration.** `mc_lower` and `mc_upper` bracket
   the finite Monte Carlo mean, conditional on accurate evidence values.
   The default uses the upper value, capped by Ville.
3. **Monte Carlo sampling.** `mcse` is a pointwise standard-error diagnostic.
   It is not a simultaneous confidence band or a correction guaranteeing
   anytime validity. Reducing quadrature tolerance does not reduce this
   source of error.

When the numerical calibration bracket straddles `alpha`, the package
reports the affected prefixes in `fit$calibration_ambiguous_n`. Reduce
`tolerance` to resolve numerical ambiguity. This does not resolve
uncertainty due to the finite null sample.

To use the deterministic population-calibration implementation instead:

```r
fit_deterministic <- ev_ttest(
  x, prior = prior, alpha = 0.05,
  calibration_method = "deterministic", tolerance = 1e-6
)
```

This option uses analytic formulas where implemented and otherwise
numerical upper/lower bounds for the null expectation. Its gap concerns
numerical integration, rather than Monte Carlo sampling. These bounds also
assume accurate floating-point evidence evaluations.

The lower-level centered two-sided helpers `ev_log_evidence()`,
`ev_calibrate()` and `ev_path()` retain their deterministic behavior.
Use `ev_ttest()` and `ev_calibration()` for the complete prior interface and
explicit calibration-method selection.

## Summaries, plots and export

`summary(fit)` reports the calibration method, final running values, first
rejection times and diagnostics. `plot(fit)` compares eVille and Ville and
marks their first observed rejections. `as.data.frame(fit)` returns the
complete trajectory.

```r
plot(fit, method = "eVille")
write.csv(summary(fit)$results, "eville_summary.csv", row.names = FALSE)
write.csv(as.data.frame(fit), "eville_trajectory.csv", row.names = FALSE)

pdf("eville_t_test.pdf", width = 7, height = 4.5)
plot(fit)
dev.off()
```

The principal trajectory columns are `n`, `log_evidence`, `p_eville`,
`p_ville`, `calibration_lower` and `calibration_gap`. For Monte Carlo results,
`calibration_lower` and `calibration_gap` describe the numerically bracketed
value after the Ville cap. The raw `mc_lower` and `mc_upper` fields bracket
the empirical mean before that cap; `mc_estimate`, `mc_gap`, `mcse` and
`mc_samples` retain its additional diagnostics. None of these empirical
brackets is a confidence interval for the population calibration.

For a proper prior, an `NA` first rejection time means no rejection in the
observed sequence. For Lebesgue measure, Ville entries are always `NA`
because that calibration is unavailable.

## Assumptions and scope

The working model is independent Gaussian observations with an unknown,
common variance. Supply observations in their sampling order. Choosing an
order after inspecting the data is not covered by the theoretical guarantee.
Every monitored prefix must have positive sample variance; missing or
nonfinite observations produce an error.

Fix the prior and test direction before observing the target data.
Parameters estimated using independent source data may be supplied through
`ev_prior()`. Robustness to prior-data conflict does not remove the Gaussian
or independence assumptions.

The package implements **one-sample t-tests with `n0 = 2`**. Two-sample and
regression tests, other monitoring start times, automatic source-prior
fitting and incremental streaming updates are not implemented. Exact
anytime validity refers to the theoretical calibration under the model;
finite Monte Carlo calibration is explicitly an approximation.

## Help, examples and citation

```r
help("eVille", package = "eVille")
?ev_ttest
?ev_prior
?ev_calibration
citation("eVille")
```

Runnable scripts are included in `inst/examples/quick-start.R` and
`inst/examples/compare-priors.R`. For example:

```r
source(system.file("examples", "compare-priors.R", package = "eVille"))
```

The notebook and scripts reproducing the paper's figures remain in the
repository's separate [`experiments/` directory](../experiments).

## Development checks

From the package source directory:

```r
devtools::document()
devtools::check()
```

Tests cover analytic reference cases, Python evidence and Monte Carlo
references, scale invariance, directional identities and agreement with
`BayesFactor::ttestBF()` under matching centered two-sided Rouder settings.
The BayesFactor comparison concerns mixture evidence, not calibrated
p-values. Null simulations do not replace the mathematical validity proof.

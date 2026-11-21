# eVille: Anytime-Valid One-Sample t-Tests

**eVille provides anytime-valid one-sample t-tests for a Gaussian mean
with unknown variance.** Specify a prior on the standardized effect,
supply observations in their sampling order, and obtain running p-values,
first rejection times, summaries, and plots.

The p-values are calibrated for sequential monitoring: the probability
of ever rejecting a true null is controlled under the stated model.
This differs from repeatedly applying an ordinary fixed-sample t-test.
Monitoring starts at the second observation.

The current implementation supports:

- Gaussian and Cauchy effect-size priors, including shifted priors.
- The centered Rouder default, Cauchy(0, 1/sqrt(2)).
- Two-sided and one-sided alternatives.
- Full and half-line Lebesgue benchmarks with their separate calibration.
- Summary tables, running p-value plots, and trajectory export.

## Installation

Install from the `r-package` directory of the GitHub repository:

```r
install.packages("remotes")  # Run once if needed.
remotes::install_github(
  "DebolinaPaul/Robust_anytime_valid_t_test",
  subdir = "r-package",
  upgrade = "never"
)
```

To install a particular archived version, also supply `ref` with its
Git tag or commit SHA.

Alternatively, from the local package directory containing `DESCRIPTION`, run:

```r
# Install devtools first if needed: install.packages("devtools")
devtools::install(".", upgrade = FALSE)
```

Then load the installed package:

```r
library(eVille)
```

## Citation and paper experiments

Use `citation("eVille")` for the software citation.

This package provides the reusable R test interface. The scripts and
notebook that reproduce the paper's simulation and DepMap figures are in
the repository's separate [`experiments/` directory](../experiments).
Those experiments also include Student-t, Laplace, and horseshoe priors;
these three families are not implemented in this R release.

The R package uses analytic formulas or deterministic numerical integration
for calibration. The simulation notebook uses Monte Carlo calibration
where specified there. Both evaluate the same theoretical calibration
function; their numerical approximations need not give identical values.

## Example

The following observations are synthetic. Their order is preserved.
The illustrative prior is specified in advance, rather than fitted to `x`.

```r
library(eVille)

x <- c(-0.3, -1.1, -0.7, -0.5, -1.4, -0.8, -0.2, -1.0)

prior <- ev_prior(
  family = "gaussian",
  location = -4,
  scale = 1,
  alternative = "less"
)

fit <- ev_ttest(x, prior = prior, mu = 0, alpha = 0.05)

summary(fit)
plot(fit)

# Numeric summary and complete results after each observation.
summary(fit)$results
results <- as.data.frame(fit)
head(results)
```

Here the null is that the population mean is nonnegative, and the
alternative is that it is negative. The prior is N(-4, 1), conditioned on
negative standardized effects. In this example, both methods first reject
at observation 8. A lower rejection threshold need not produce an earlier
crossing on every sequence.

## Summaries, plots, and export

`summary(fit)` reports the final running p-values, first rejection times,
and calibration diagnostics. `plot(fit)` compares eVille and Ville and
marks their first observed rejections.

```r
plot(fit, method = "eVille")
plot(fit, main = "Anytime-valid one-sample t-test")

# These commands create files in your current working directory.
write.csv(summary(fit)$results, "eville_summary.csv", row.names = FALSE)
write.csv(as.data.frame(fit), "eville_trajectory.csv", row.names = FALSE)

pdf("eville_t_test.pdf", width = 7, height = 4.5)
plot(fit, main = "Anytime-valid one-sample t-test")
dev.off()
```

## Choosing a prior and alternative

The prior is on the standardized effect `theta = (mu - mu0) / sigma`.
The `mu` argument of `ev_ttest()` specifies the null boundary `mu0`.

```r
# Shifted two-sided Gaussian prior
p_gaussian <- ev_prior("gaussian", location = 3, scale = 1)

# Shifted two-sided Cauchy prior
p_cauchy <- ev_prior("cauchy", location = 3, scale = 1)

# Fixed centered Rouder prior
p_rouder <- ev_prior("rouder")

# Lebesgue measure on negative effects
p_flat <- ev_prior("lebesgue", alternative = "less")
```

| `alternative` in `ev_prior()` | Null hypothesis | Alternative |
|---|---|---|
| `"two.sided"` | mu = mu0 | mu != mu0 |
| `"less"` | mu >= mu0 | mu < mu0 |
| `"greater"` | mu <= mu0 | mu > mu0 |

For proper one-sided priors, the specified distribution is conditioned on
the selected half-line. Location and scale describe the distribution before
conditioning. Gaussian scale is its standard deviation; Cauchy scale is not
a standard deviation.

Fix the prior independently of the target observations. Priors estimated
from independent source data can be supplied through their fitted parameters;
automatic source-prior fitting is not currently implemented.

## Results and the Lebesgue benchmark

`as.data.frame(fit)` contains:

| Column | Meaning |
|---|---|
| `n` | Number of observations used |
| `log_evidence` | Natural log of the current mixture evidence |
| `p_eville` | Running eVille value used for rejection |
| `p_ville` | Running ordinary-Ville value |
| `calibration_lower` | Lower calibration bound at the running evidence |
| `calibration_gap` | Difference between the upper and lower calibration values |

For proper priors, `NA` in `fit$stopping_times` means no rejection within the
observed sequence. Ordinary Ville is unavailable for Lebesgue measure, so its
Ville entries are always `NA`.

```r
flat <- ev_ttest(x, prior = p_flat, alpha = 0.05)
flat$stopping_times
```

## Numerical precision

General shifted-prior calibration uses an upper approximation to the ideal
extended Ville calibration, together with a reported integration gap. The
upper bound is conservative for exact evidence evaluations. Centered
two-sided Gaussian and full/half-Lebesgue calibrations use analytic formulas.

If the calibration bracket straddles `alpha`, the test warns and records the
affected prefixes in `fit$calibration_ambiguous_n`. Reduce the tolerance to
resolve the decision more precisely:

```r
precise_fit <- ev_ttest(
  x, prior = prior, alpha = 0.05, tolerance = 1e-6
)
```

The default tolerance is `1e-4`. The reported gap concerns null-angle
integration; it is not a rigorous bound on every floating-point operation.
Evidence quadrature is checked by refinement and an omitted-tail bound.
Numerical tests do not replace the mathematical validity argument.

## Assumptions and current scope

The working model is

$$
Y_1,Y_2,\ldots \stackrel{\mathrm{iid}}{\sim}N(\mu,\sigma^2),
\qquad \sigma>0\text{ unknown}.
$$

Supply data in their original sampling order. Sorting the observed values
or selecting the most favourable order is not covered by the guarantee.
Robustness to prior-data conflict does not remove the Gaussian or
independence assumptions. Each monitored prefix must have positive sample variance.
Missing or nonfinite observations produce an error and are not silently
removed. The prior and test direction are fixed before monitoring.

The package currently implements **one-sample t-tests with `n0 = 2`**.
Two-sample and regression tests, other monitoring start times, additional
prior families, automatic source-prior fitting, and incremental streaming
updates are not yet implemented.

## Help

```r
help("eVille", package = "eVille")
?ev_ttest
?ev_prior
```

A runnable example is included in `inst/examples/quick-start.R` in the source.
After installation, it can be run with:

```r
source(system.file("examples", "quick-start.R", package = "eVille"))
```

## Development checks

From the package source directory:

```r
devtools::document()
devtools::check()
```

The tests include analytic reference cases, Python comparisons, scale
invariance, directional identities, calibration bounds, and agreement with
`BayesFactor::ttestBF()` under matching centered two-sided Rouder settings.
That comparison concerns mixture evidence, not calibrated p-values.

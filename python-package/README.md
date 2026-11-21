# eVille: anytime-valid one-sample t-tests in Python

Python implementation accompanying **"Bayes-assisted anytime-valid t-tests
for Gaussian means with unknown variance"**, by Debolina Paul, Patrick
Rebeschini and François Caron.

**eVille implements one-sample t-tests for a Gaussian mean with unknown
variance**, using Ville and extended Ville calibration. Supply observations
in their sampling order to obtain running values, first rejection times,
summaries and plots. Monitoring starts at observation 2.

The package supports Gaussian, Cauchy, Student-t5, Laplace and genuine
horseshoe effect-size priors, together with the centered Rouder prior and
full or half-line Lebesgue benchmarks. It uses the prior parameterizations
and shared Monte Carlo null sample from the accompanying simulation notebook.

**Calibration:** the default for proper priors is a finite Monte Carlo
approximation. The exact anytime-valid guarantee belongs to the population
calibration; the empirical approximation does not itself establish that
guarantee. Numerical and Monte Carlo diagnostics are reported separately.
Deterministic population calibration is also available.

## Installation

Python 3.10 or later is required. From the repository root, install the
local package:

```bash
python -m pip install ./python-package
```

From GitHub, it can also be installed with:

```bash
python -m pip install "git+https://github.com/DebolinaPaul/Robust_anytime_valid_t_test.git#subdirectory=python-package"
```

The import name is `eville`. This installation uses GitHub; it does not
require a PyPI release. See [the setup guide](docs/setup.md) for a Mac virtual
environment and development checks, and the [interface guide](docs/api.md)
for functions and result methods.

## A one-sample t-test

The observations below are synthetic. The prior is specified independently
of these observations.

```python
from eville import ev_prior, ev_ttest
import matplotlib.pyplot as plt

x = [-0.3, -1.1, -0.7, -0.5, -1.4, -0.8, -0.2, -1.0]
prior = ev_prior("gaussian", location=-4, scale=1, alternative="less")
fit = ev_ttest(x, prior=prior, mu=0, alpha=0.05)

print(fit.summary())
print(fit.to_dataframe().head())
ax = fit.plot()
plt.show()
```

This tests a nonnegative mean against a negative mean. The prior is
Normal(-4, 1), conditioned on negative standardized effects. Ville and
eVille use the same mixture evidence; their calibration differs.

Export results and figures using the returned pandas and Matplotlib objects:

```python
fit.summary().to_csv("eville_summary.csv", index=False)
fit.to_dataframe().to_csv("eville_trajectory.csv", index=False)
ax.figure.savefig("eville_t_test.pdf", bbox_inches="tight")
```

For your own CSV, select the numeric outcome column in its observation order:

```python
import pandas as pd

x = pd.read_csv("observations.csv")["response"].to_numpy()
fit = ev_ttest(x, prior=prior, mu=0, alpha=0.05)
```

Replace `response` with your column name and specify the prior independently
of these target observations. Missing values cause an error rather than being
silently removed.

## Priors and hypotheses

The prior is on `theta = (population_mean - mu) / sigma`, where `mu` is the
null mean passed to `ev_ttest()`. Write `m = location` and `s = scale`.

| Family | Distribution before one-sided conditioning | Default scale |
|---|---|---:|
| `"gaussian"` | Normal with mean `m` and standard deviation `s` | 1 |
| `"cauchy"` | Cauchy with location `m` and scale `s` | `1 / sqrt(2)` |
| `"student_t5"` | `m + s * T`, with `T` distributed as Student-t with 5 degrees of freedom | 1 |
| `"laplace"` | Density `exp(-abs(theta - m) / s) / (2 * s)` | 1 |
| `"horseshoe"` | `m + s * lambda * Z`, with independent standard normal `Z` and half-Cauchy(0, 1) `lambda` | 1 |

The Student-t5 variance is `5 * s**2 / 3`; the Laplace variance is
`2 * s**2`. The horseshoe has an infinite density at `m`, and its scale is
not a standard deviation. It is the normal--half-Cauchy mixture above,
not a logarithmic approximation to that density.

```python
from eville import ev_prior_density

prior = ev_prior("student_t5", location=3, scale=1)
print(ev_prior_density([-1, 0, 1], prior))

rouder = ev_prior("rouder")
flat = ev_prior("lebesgue", alternative="less")
```

`"rouder"` is a centered Cauchy with default scale `1 / sqrt(2)`.
`"lebesgue"` is an improper mixing measure. Its full or half-line evidence
has a separate analytic extended Ville calibration; ordinary Ville is
unavailable for this benchmark.

| `alternative` | Null hypothesis | Alternative hypothesis |
|---|---|---|
| `"two.sided"` | mean = mu | mean != mu |
| `"less"` | mean >= mu | mean < mu |
| `"greater"` | mean <= mu | mean > mu |

For a one-sided test, each proper prior is conditioned on the corresponding
half-line. Its location and scale describe the distribution before
conditioning. The notebook uses two-sided priors; one-sided support is an
additional package feature.

## Calibration and agreement with the notebook

Proper priors use `calibration_method="monte_carlo"` by default. The
package bundles the same 100,000 null t statistics derived from independent
standard-normal pairs in the notebook, using NumPy seed `20260510`.
Loading this fixed sample does not alter the user's random-number generator.

The empirical calibration averages `min(1, M2_null / c)` over those draws.
Adaptive numerical bounds evaluate this finite average. As in the notebook,
the upper numerical value is capped by Ville's bound `min(1, 1/c)`.
Thus, eVille's computed running values are no larger than Ville's for the
same proper mixture. This cap does not remove Monte Carlo sampling error.

The default absolute calibration tolerance is `1e-5`. Trajectories also
use relative tolerance `1e-3`, matching the notebook's paper configuration.
Matching the notebook or R package requires the same observations, prior,
alternative, null draws and tolerances; floating-point computations can
still differ within numerical accuracy.

Calibration objects cache calculations and can be reused with the same
prior:

```python
from eville import ev_calibration, ev_null_sample

null_t = ev_null_sample()
prior = ev_prior("student_t5", location=3)
calibration = ev_calibration(prior, null_t=null_t)
fit = ev_ttest(x, prior=prior, calibration=calibration)
```

A custom `null_t` sample must contain independent draws from the null
Student-t distribution with one degree of freedom, generated independently
of the target data. It changes the Monte Carlo approximation. For another
prior, construct a new calibrator; the same null sample can be reused.

To calculate the population calibration without Monte Carlo sampling:

```python
fit_deterministic = ev_ttest(
    x,
    prior=prior,
    calibration_method="deterministic",
    tolerance=1e-6,
)
```

This method uses analytic formulas or numerical integration, with numerical
diagnostics. Lebesgue calibration is analytic regardless of the selected
proper-prior method. See [numerical methods](docs/numerics.md) for the
distinction between quadrature uncertainty and Monte Carlo error.

## Assumptions and scope

- The model assumes independent Gaussian observations with an unknown,
  common variance. Robustness to prior-data conflict does not remove
  these assumptions.
- Supply the observations in their sampling order. Choosing their order
  after inspecting their values is not covered by the theoretical guarantee.
- Specify the prior and direction before observing the target data.
  Parameters learned from independent source data can be passed to `ev_prior()`.
- Monitoring starts at `n0=2`. Every monitored prefix must have positive
  sample variance. Missing and nonfinite values produce an error; they are
  not silently discarded.
- This is a **one-sample t-test** package. For paired data, applying it to
  independent Gaussian within-pair differences is a one-sample analysis.
  Independent two-sample tests, regression tests, other start times,
  automatic prior fitting and incremental updates are not implemented.

The methods return the full trajectory even if rejection occurs before the
last observation. A missing first rejection time means no rejection within
the supplied observations; an unavailable Ville result for Lebesgue is a
different case and is labelled accordingly.

## Examples and development checks

From `python-package/` after installation:

```bash
python examples/quick_start.py
python examples/compare_priors.py
python -m pip install -e ".[dev]"
python -m pytest
python -m build
```

Examples use synthetic observations and save results in `example_outputs/`.
Tests check the numerical engine, analytic cases, notebook references,
invariance and result methods. Such checks diagnose implementation errors;
they do not replace the mathematical validity argument.

The main-paper and supplementary simulation figures remain in
[`experiments/simulations/`](../experiments/simulations/), and the DepMap
analysis remains in [`experiments/depmap/`](../experiments/depmap/).

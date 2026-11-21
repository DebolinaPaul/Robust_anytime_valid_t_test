# Python interface

Import public functions from `eville`. Internal modules beginning with an
underscore contain implementation details and are not the public interface.

| Function or method | Purpose |
|---|---|
| `ev_prior(...)` | Specify a fixed effect-size prior or Lebesgue measure |
| `ev_prior_density(theta, prior, log=False)` | Evaluate its density or mixing intensity |
| `t_statistic(x, mu=0)` | Calculate the ordinary sample t statistic |
| `ev_log_evidence(x, prior=None, mu=0)` | Calculate log evidence for one prefix |
| `ev_log_mixture(t, n, prior)` | Calculate log evidence from a t statistic and sample size |
| `ev_null_sample(n=100000)` | Load the shared null t statistics |
| `ev_calibration(prior, ...)` | Construct and cache a calibration for one prior |
| `ev_pvalues(...)` | Evaluate calibration at supplied log-evidence thresholds |
| `ev_ttest(x, ...)` | Calculate the complete ordered one-sample trajectory |
| `fit.summary()` | Return a pandas table of decisions and first rejection times |
| `fit.to_dataframe()` | Return a copy of the trajectory and diagnostics |
| `fit.plot(...)` | Plot stored running values and return a Matplotlib Axes |

Use Python's `help()` for signatures and details:

```python
from eville import ev_prior, ev_ttest, ev_calibration, ev_pvalues

help(ev_prior)
help(ev_ttest)
help(ev_calibration)
help(ev_pvalues)
```

## The main test

```python
fit = ev_ttest(
    x,
    prior=None,
    mu=0.0,
    alpha=0.05,
    n0=2,
    calibration=None,
    tolerance=1e-5,
    calibration_method="monte_carlo",
    null_t=None,
    relative_tolerance=1e-3,
)
```

- `x` is a finite numeric one-dimensional sequence in observation order.
- `prior=None` selects the default centered Gaussian prior with scale 1.
  Use `ev_prior(...)` to choose a different family, location, scale or direction.
- `mu` is the null mean, rather than a prior location. The prior is on the
  standardized difference from this mean.
- `alpha` is strictly between zero and one. `n0` must be 2.
- `calibration_method` selects `"monte_carlo"` or `"deterministic"` for
  proper priors. Lebesgue uses its analytic formula.
- `tolerance` controls numerical calibration accuracy. The trajectory
  additionally uses `relative_tolerance`; `None` requests absolute
  tolerance only. Neither controls finite Monte Carlo sampling error.
- A supplied `calibration` must match the prior. Its method and tolerance
  take precedence over these settings; do not also supply `null_t`.

The returned result stores `prior`, `alpha`, `mu`, `calibration`,
`stopping_times` and `calibration_ambiguous_n`. The last field identifies
prefixes where the numerical calibration bracket straddles `alpha`.
Tighter numerical tolerances can resolve that ambiguity, but do not remove
Monte Carlo sampling uncertainty.

`fit.stopping_times` maps `"eVille"` and `"Ville"` to their first observed
rejection times, or `None`. Since `None` also represents unavailable Ville
calibration for Lebesgue, use the summary's `decision` column to distinguish
the two cases.

## Result tables

`fit.summary()` contains `method`, `final_running_p`, `first_rejection`,
`decision`, `calibration`, `mc_samples` and `final_mcse`. Additional settings
are available in its pandas `attrs` dictionary. CSV export does not retain
those attributes, so retain the analysis script alongside exported tables.

`fit.to_dataframe()` contains evidence, running values and calibration
diagnostics at each prefix. The Monte Carlo diagnostics refer to the
empirical calibration, not confidence intervals for the population
calibration; see [numerical methods](numerics.md).

## Plots

```python
ax = fit.plot(
    method="both",       # "both", "eVille" or "Ville"
    ax=None,             # Or an existing Matplotlib Axes
    log_y=False,
    mark_rejections=True,
    title=None,
    legend=True,
)
```

Plots use stored results and do not rerun the test. Lebesgue plots contain
only eVille; requesting Ville alone for Lebesgue is an error. Monte Carlo
eVille is labelled explicitly. For a logarithmic vertical axis, floating-point
zeros are moved to a positive display floor without changing stored values.

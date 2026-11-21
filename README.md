# Bayes-assisted anytime-valid t-tests for Gaussian means with unknown variance

**Debolina Paul, Patrick Rebeschini and François Caron**

This repository contains statistical software and numerical experiments
for anytime-valid one-sample testing of a Gaussian mean with unknown
variance. The methods incorporate prior information about a standardized
effect size and compare reciprocal Ville calibration with extended Ville
calibration. The principal focus is the effect of prior-data conflict on
power and stopping times.

R and Python implementations are provided, together with the simulation
studies and DepMap analyses.

## Statistical framework

For independent observations

\[
Y_i \sim N(\mu,\sigma^2), \qquad \sigma>0,
\]
the standardized effect is
\[
\theta=\frac{\mu-\mu_0}{\sigma}.
\]

The two-sided test considers \(H_0:\mu=\mu_0\) against
\(H_1:\mu\ne\mu_0\). One-sided alternatives are also implemented.
Monitoring begins at the second observation.

For a fixed proper prior, Ville and extended Ville use the same mixture
evidence process and differ in its calibration. The available prior
families are Gaussian, Cauchy, Student-t with five degrees of freedom,
Laplace and horseshoe. The implementations also include the centered
Rouder Cauchy prior with scale \(1/\sqrt{2}\), and full-line or half-line
Lebesgue benchmarks.

Priors must be specified independently of the target observations.
They may be estimated from independent source data. Robustness to
prior-data conflict does not remove the Gaussian sampling and
independence assumptions.

## Calibration

For proper priors, the R and Python packages use Monte Carlo calibration
by default, with the same 100,000 simulated null statistics used in the
simulation notebook. Numerical integration bounds and Monte Carlo
standard errors are reported separately.

The finite Monte Carlo calculation approximates the population
calibration. The exact anytime-valid guarantee concerns the population
calibration; the empirical numerical bounds do not bound Monte Carlo
sampling error.

Deterministic population calibration is also available in both packages.
The improper Lebesgue benchmarks use separate analytic extended Ville
calibrations. Ordinary Ville calibration is unavailable for these
improper measures.

## Repository structure

| Directory | Contents |
|---|---|
| [`r-package/`](r-package/) | R package source, documentation, examples and tests |
| [`python-package/`](python-package/) | Python package source, documentation, examples and tests |
| [`experiments/simulations/`](experiments/simulations/) | Simulation notebook, saved numerical results and figures |
| [`experiments/depmap/RPL38/`](experiments/depmap/RPL38/) | RPL38 analysis with the recorded Female sex restriction |
| [`experiments/depmap/RPS23/`](experiments/depmap/RPS23/) | RPS23 analysis with the recorded Male sex restriction |
| [`data/`](data/) | DepMap metadata, data provenance and download instructions |

Experiment outputs are stored within their respective experiment
directories.

## R package

Install from GitHub:

```r
install.packages("remotes")

remotes::install_github(
  "DebolinaPaul/Robust_anytime_valid_t_test",
  subdir = "r-package",
  upgrade = "never"
)

library(eVille)
```

The principal interface is `ev_ttest()`, with priors specified through
`ev_prior()`. Results support `summary()`, `plot()` and
`as.data.frame()`.

See the [R package documentation](r-package/README.md) for prior
parameterizations, calibration options and examples.

## Python package

Python 3.10 or later is required. Install in an appropriate Python
environment:

```bash
python -m pip install --upgrade pip

python -m pip install \
  "git+https://github.com/DebolinaPaul/Robust_anytime_valid_t_test.git#subdirectory=python-package"
```

The import name is `eville`. The principal interfaces are `ev_prior()`
and `ev_ttest()`. Results provide `summary()`, `plot()` and
`to_dataframe()` methods.

See the [Python package documentation](python-package/README.md) for
examples and numerical details.

## Simulation studies

The simulations examine sensitivity to prior location under Gaussian
sampling with mean 0.5 and variance 1.

- Fixed-time rejection probabilities use 10,000 simulated streams,
  sample sizes 20, 50, 100 and 150, and significance level 0.1.
- Running calibrated values use one independent trajectory of
  150 observations and significance level 0.05.
- The main figures include Gaussian, Cauchy and Student-t5 priors.
  Supplementary figures additionally include Laplace and horseshoe priors.

Fixed-time rejection probability and rejection by a monitoring horizon
are distinct quantities; the power figures report the former.

To reproduce the simulations, run the following commands from the
repository root:

```bash
python -m pip install numpy scipy pandas matplotlib jupyterlab
cd experiments/simulations
python -m jupyter lab Code_simulation.ipynb
```

Set `RUN_MODE = "paper"` and select **Restart Kernel and Run All**.
The notebook is self-contained and does not require either package.

Results are written to `simulation_outputs/paper_<configuration hash>/`.
Each run records its configuration, software versions, simulated
observations, calibration diagnostics and numerical results. The
principal figures are:

- `power_vs_m_priors.pdf`
- `pvalue_vs_n_simulated.pdf`
- `power_vs_m_priors_supplement.pdf`
- `pvalue_vs_n_simulated_supplement.pdf`

## DepMap analyses

The real-data analyses use CRISPR gene-knockout measurements from
[DepMap Public 24Q4](https://doi.org/10.25452/figshare.plus.27993248.v1).
Negative Chronos gene-effect scores indicate reduced cell growth
following knockout.

Each eligible cancer lineage is held out in turn. Priors are fitted
using the remaining eligible lineages for the same gene. The analysis
restricts recorded sex, retains one model per recorded patient and
requires at least ten observations in each source and target lineage.

Each target uses 100 random observation orders, shared across methods.
An ascending ModelID order is retained separately as an illustration.
Permutation summaries describe sensitivity to observation order
conditional on the observed dataset; they are not estimates of
repeated-sample power.

The experiment directories include single-gene extracts sufficient
for their respective analyses. Instructions for downloading the full
gene-effect matrix and analysing other genes are provided in
[`data/README.md`](data/README.md).

The following commands reproduce the significance levels recorded in
the currently committed results:

```bash
python -m pip install -r experiments/depmap/RPL38/requirements.txt

python -u experiments/depmap/RPL38/run_analysis.py \
  --alpha 0.005 \
  --min-source-n 10 \
  --min-target-n 10 \
  --permutations 100 \
  --out experiments/depmap/RPL38/results_reproduced

python -u experiments/depmap/RPS23/run_analysis.py \
  --alpha 0.001 \
  --min-source-n 10 \
  --min-target-n 10 \
  --permutations 100 \
  --out experiments/depmap/RPS23/results_reproduced
```

Each completed run records its parameters in `settings.json`.
Changing the significance level or cohort restrictions defines a
different analysis and requires recomputing the results.

Saved outputs include fitted priors, stopping times, summaries,
illustrative trajectories, figures and validation reports. Saved
results can be checked independently:

```bash
python experiments/depmap/RPL38/validate.py \
  --results experiments/depmap/RPL38/results

python experiments/depmap/RPS23/validate.py \
  --results experiments/depmap/RPS23/results
```

The gene and illustrative target selections are exploratory.
The analyses do not provide multiplicity or selective-inference
adjustments for those selections.

## Software checks

From a local repository clone, run the Python tests with:

```bash
python -m pip install -e "./python-package[dev]"
python -m pytest python-package/tests
```

Run the R package checks with:

```r
devtools::check("r-package")
```

The tests compare numerical calculations with analytic formulas,
independent integrations and cross-language reference values.
They also check scale invariance, directional identities and agreement
with `BayesFactor::ttestBF()` under matching two-sided Rouder settings.

Numerical checks assess the implementation within the tested cases.
They do not replace the mathematical validity argument or establish
that a real dataset satisfies the sampling assumptions.

## Citation

When using these methods, cite:

Paul, D., Rebeschini, P. and Caron, F.
*Bayes-assisted anytime-valid t-tests for Gaussian means with unknown variance.*

For analyses using DepMap, also cite:

Broad DepMap (2024). *DepMap 24Q4 Public*.
https://doi.org/10.25452/figshare.plus.27993248.v1.

## Licence

The software is distributed under the [MIT licence](LICENSE).
DepMap data remain subject to their source licence and attribution
requirements, as described in [`data/README.md`](data/README.md).

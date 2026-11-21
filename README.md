# Robust Anytime-Valid t-Tests

Code accompanying **“Bayes-assisted anytime-valid t-tests for Gaussian
means with unknown variance”**, by Debolina Paul, Patrick Rebeschini,
and François Caron.

This project studies anytime-valid one-sample t-tests that incorporate
prior information about a standardized effect size. We compare ordinary
Ville calibration with extended Ville calibration, focusing on the
effect of prior-data conflict on power and stopping times.


## Statistical setting

We consider independent Gaussian observations with unknown variance,

$$
Y_i \stackrel{\mathrm{iid}}{\sim} N(\mu,\sigma^2),
\qquad \sigma>0,
$$

and place a working prior on the standardized effect

$$
\theta=\frac{\mu-\mu_0}{\sigma}.
$$

The main theoretical analysis concerns testing
$H_0:\mu=\mu_0$ against $H_1:\mu\ne\mu_0$.
The implementation also supports one-sided alternatives.

The tests produce p-values calibrated for sequential monitoring.
Under the stated model, the probability of ever rejecting a true null
is controlled at the chosen significance level.

The prior is fixed before examining the target observations. It can
incorporate information from independent source data. Robustness to
prior-data conflict does not remove the Gaussian or independence
assumptions.

## R package: eVille

**eVille: Anytime-Valid One-Sample t-Tests** provides:

- Gaussian and Cauchy effect-size priors, including shifted priors.
- The centered Rouder default, Cauchy with scale $1/\sqrt{2}$.
- One-sided and two-sided alternatives.
- Ville and extended Ville calibration for proper prior mixtures.
- Separately calibrated full and half-line Lebesgue benchmarks.
- Running p-values, first rejection times, summary tables, and plots.

The current implementation begins monitoring at the second observation.
Ordinary Ville calibration is unavailable for the improper Lebesgue
benchmark.

The package source will be provided in `eVille/`. GitHub installation
instructions will be added once it has been uploaded.

## Publication code organisation

The intended organisation is:

| Location | Purpose |
|---|---|
| `eVille/` | R package source, documentation, and tests |
| `experiments/simulations/` | Simulation code, configurations, and plotting scripts |
| `experiments/depmap/` | DepMap preprocessing and real-data analysis |
| `data/README.md` | Data sources and download instructions |
| `results/` | Numerical results and figures used in the paper |
| `reproduce/` | Reproduction instructions and figure/table manifest |
| `JRSSB/` | Main manuscript and supplementary material |
| `archive/` | Historical experiments and earlier versions |

The Python experiment code will include the additional prior families
used in the paper, beyond those currently supported by the R package.

## Experiments

### Simulations

The simulations examine how prior location affects rejection
probabilities and running p-values, comparing Ville, extended Ville,
and the Lebesgue benchmark.

Fixed-time rejection probabilities and first rejection times are
different summaries and are reported separately.

### DepMap application

The real-data experiment uses RPL38 CRISPR gene-knockout responses
from DepMap Public 24Q4. Cancer lineages are held out in turn, and
prior information is constructed using the remaining source lineages.

The required input files are:

- `CRISPRGeneEffect.csv`
- `Model.csv`

Original data release:
[DepMap 24Q4 Public](https://doi.org/10.25452/figshare.plus.27993248.v1).

The analysis specifies the filtering,
patient deduplication, prior fitting, and observation-order procedures
needed to reproduce the reported results.

## Reproducibility

The publication release will include:

- Exact simulation and analysis configurations.
- Random seeds and observation-order information.
- Software dependencies and version information.
- Commands to reproduce the main and supplementary results.
- A manifest linking each figure and table to its generating script.
- Expected outputs and approximate runtimes.

Numerical results will be saved separately from plotting, allowing
figures to be regenerated without repeating the full calculations.

## Citation

A versioned software citation and an archived-release DOI will be
provided with the publication release.

## Licence

Project software is distributed under the [MIT licence](LICENSE).
Third-party templates retain their original licence notices.
Source datasets remain subject to their original terms.

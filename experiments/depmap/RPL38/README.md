# RPL38: one final source-to-target analysis

This package contains one implementation and one fitting analysis. It uses the uploaded DepMap Public 24Q4 RPL38 data and the recorded Female sex filter. It does not contain an environment, old scripts, an IQR comparison run, or a directory for every held-out lineage.

## Run on macOS / VS Code

Use Python 3.10 or newer. In the VS Code terminal, change into this extracted folder, then run:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -u run_analysis.py --permutations 100 --overwrite
```

The ZIP already includes a completed `results/` folder. `--overwrite` replaces only a recognized completed output folder created by this package; calculations and validation finish in a temporary folder before replacement. To retain the bundled results, use `--out results_new` instead. A failed run leaves its temporary folder for diagnosis and does not replace an existing completed result.

To redraw figures without refitting, run:

```bash
python figures.py --results results
```

To independently check the saved calculations, run:

```bash
python validate.py --results results
```

## What is tested

For each lineage, Y is the raw Chronos score. The working model is iid N(mu, sigma^2), with H0: mu >= 0 versus H1: mu < 0. Negative scores indicate reduced growth after knockout. Computation uses X=-Y and a positive standardized alternative. Monitoring begins at n0=2, with alpha=0.05. This package implements calibration only for n0=2.

Models must have finite scores and complete lineage and PatientID metadata. Recorded Sex is filtered first, then the lexicographically first ModelID is retained for each PatientID. There are 451 retained donors, 26 lineages, and 24 testable lineages with at least two observations. The 17 lineages with N>=10 form the main cohort; the other seven are flagged as small-group sensitivity results. The two singleton lineages are listed as untestable.

Every target lineage is excluded in full from its source data. Sources need at least ten retained observations. One standardized sample mean per eligible source lineage enters the fit, with equal lineage weight. Source-estimate uncertainty is not modelled. Source and target PatientID sets are checked to be disjoint, which does not establish biological independence.

## One fitting analysis

Gaussian, Cauchy, Laplace and Student-t5 full-family location/scale parameters are fitted by maximum likelihood, then the fitted density is conditioned on the alternative half-line. This is not a truncated-family MLE. Cauchy and t5 fits use multistart numerical optimization and are not certified global maxima.

**Horseshoe exception:** its unrestricted location likelihood is unbounded, so a finite free-location MLE does not exist. It retains the explicitly labelled median/IQR fit and is shown as Horseshoe*. There is no second IQR analysis. The IQR also supplies optimizer initial values; this does not turn the other four final fits into IQR fits.

Rouder is the fixed Cauchy(0,1/sqrt(2)) conditioned on the alternative half-line. The flat comparator is genuinely improper Lebesgue measure, evaluated using extended Ville only. All of its ordinary Ville and paired-comparison fields are blank.

## Repeated orders and summaries

Each lineage uses the fixed ascending ModelID order 0 plus 100 random orders. Random orders are identical across priors and calibration methods, and stable per-lineage seeds make results independent of lineage execution order. Summaries exclude order 0.

Nonrejection by N is recorded as missing, not as a rejection at N+1. Restricted stopping means average min(tau,N), with tau=infinity for nonrejection. The empirical unconditional median is the lower 0.5 quantile and is NR unless at least half the random orders reject. Both rejection fractions and restricted means remain in the table.

The paired-savings figure distinguishes: equal finite stopping times, 1/2/at least 3 observations saved when both reject, eVille-only rejection, and neither rejection. This avoids a misleading zero-saving bar for N=2 groups where eVille rejects at N but Ville does not. The flat comparator is absent from paired-savings panels.

## What the real data can establish

Random permutations of one observed sample describe sensitivity to observation order. They are not independent biological samples and do not estimate repeated-sample power. The target standardized mean is estimated, not known.

The conflict display is descriptive: abs(target estimated standardized mean - full-family prior location) / prior scale. The scale is not necessarily a standard deviation. Compare patterns within each prior panel; this is not a calibrated conflict test. Sample size, target effect and source prior vary across lineages. Therefore the scatter cannot establish a causal or monotone conflict-power relation or an asymptotic return to Lebesgue.

The real-data example supports modest eVille-versus-Ville savings. Myeloid shows loss of the benefit of fitted informative priors relative to the flat baseline, while eVille remains closer to that baseline than Ville. Keep the simulated conflict-power/limiting-baseline claim in the simulation section. A flat limit depends on the chosen sequence of priors; increasing conflict alone does not universally imply it.

The gene, sex filter and illustrative target were selected retrospectively. These plots are exploratory. Source-only fitting does not remove that selection, and the package does not apply multiple-testing adjustment across lineages or the prior candidate search.

## Compact outputs

- `summary.csv`: one row per lineage/prior; first file to read. Contains stopping means, lower medians, rejection fractions, paired outcomes, and differences from the flat baseline.
- `cohort.csv`: every retained lineage, eligibility, size and fixed-sample Gaussian effect interval.
- `prior_fits.csv`: fits, fit diagnostics and calibration thresholds for all targets.
- `stopping_times.csv.gz`: one compressed table of all fixed/random-order first rejections; no per-lineage copies.
- `critical_boundaries.csv.gz`: consolidated boundaries for independent verification; a trailing missing boundary means unnecessary after all evaluated orders rejected.
- `observations.csv.gz`: the retained, deduplicated data. Seeds and ModelID sorting reconstruct all permutation indices.
- `example_paths.csv`: full original-order p-value paths for the declared illustrative target only.
- `settings.json`, `data_filter_audit.json`, `boundary_checks.json`, `validation.json`, `manifest.json`: reproducibility and numerical checks.
- `figures/`: four final figures; no redundant all-100%-rejection heatmap.

The main and small-group results are kept in these same tables with an explicit `cohort` field. There are no separate versions of the analysis.

## Source files

`run_analysis.py` is the single analysis entry point. `cohort.py` handles selection and intervals; `fitting.py` fits priors; `sequential.py` handles paired orders and crossings; `engine.py` and `improper_flat.py` evaluate evidence/calibration; `effect_intervals.py` computes fixed-sample intervals; `figures.py` renders the consolidated results; `validate.py` independently verifies saved calculations.

`data/CRISPRGeneEffect.csv` is the RPL38-only extraction (1178 rows), not the complete approximately 429 MB DepMap matrix. `data/Model.csv` is the supplied metadata file. `data/provenance.json` records the original DepMap release DOI/download references and checksums. Input SHA256 digests are also recorded in each run.

## Numerical scope

Stable offset-Welford variance, log-domain likelihood mixtures, direct clipped-expectation calibration, the transformed horseshoe survival integral and the analytic flat comparator are used. The legacy grid calibration is removed. Fresh data analysis is compared with independent boundary and likelihood calculations; saved results are checked at higher quadrature order for the illustrative target.

These numerical checks are not interval-arithmetic error certificates or evidence that Gaussianity/independence hold. Extreme priors conditioned on a vanishingly small half-line mass can still be poorly resolved by the proper-prior quadrature; the current RPL38 fits are outside that regime. Do not extrapolate the audit to arbitrary new prior parameters.

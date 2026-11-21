# RPS23: a second DepMap example

This package repeats the uploaded RPL38 experiment for **RPS23**, restricting
the data to cell lines with recorded sex **Male**. The illustrative held-out
lineage is **Prostate**. Every eligible lineage is held out in turn; eligibility
requires at least **10** observations in both source and target lineages.

The completed analysis, input data, six PDF figures, and numerical checks are
included. There is one code version. The numerical fitting, evidence,
calibration, stopping-time, cohort, and validation modules are unchanged from
the supplied `RPL38_final.zip`.

## Run on macOS in VS Code

Open this extracted `RPS23_final` folder in VS Code, then open its terminal.
Use Python 3.10 or newer:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -u run_analysis.py --permutations 100 --min-target-n 10 --overwrite
```

`--overwrite` replaces this package's included completed results after the new
run passes validation. An unsuccessful run is retained separately for diagnosis.
To keep the included results, instead use:

```bash
python -u run_analysis.py --permutations 100 --min-target-n 10 --out results_rerun
```

The two input CSV files are already in `data/`; paths are relative to the
script, so no path edits are necessary. `data/CRISPRGeneEffect.csv` contains the
RPS23 column and all original model rows, not the full multi-gene matrix.

To recreate figures from existing results without refitting:

```bash
python figures.py --results results
```

## Results to open

- `results/figures/raw_stopping_times_all_lineages.pdf`: the same five-panel
  plot of raw mean stopping times, with fitted-prior eVille and Ville, Lebesgue
  eVille, and Rouder eVille.
- `results/figures/paired_savings_all_lineages.pdf`: paired stopping outcomes
  across eligible held-out lineages.
- `results/figures/conflict_and_flat_baseline.pdf`: stopping times and prior
  discrepancy across all lineages.
- `results/figures/conflict_and_flat_baseline_example.pdf`: the same plot
  restricted to Prostate. This replaces the old Myeloid-specific filename.
- `results/figures/example_forest_and_evidence.pdf`: Prostate forest plot and
  the Gaussian evidence trajectory in ascending ModelID order.
- `results/figures/example_running_pvalues.pdf`: running anytime-valid
  p-values in the same fixed order, for all priors.
- `results/paper_stopping_times.csv`: the requested four-column table, with
  equal-lineage averages and Prostate alone, rounded to two decimals.
- `results/summary.csv`: unrounded lineage/prior summaries.
- `results/stopping_times.csv.gz`: every first rejection time. Order 0 is the
  fixed ascending ModelID order; orders 1–100 are the random permutations used
  in the reported means. Pandas reads this compressed CSV directly.
- `results/prior_fits.csv`: fitted parameters, thresholds, and fitting checks.
- `results/cohort.csv`: eligible and excluded lineages, sizes and effect estimates.
- `results/validation.json` and `results/boundary_checks.json`: numerical checks.
- `RESULTS.md`: concise results and interpretation.

## Same experiment and methods

After filtering recorded sex, finite scores, cancer lineages and complete
patient/lineage metadata, retain the first ModelID for each PatientID in
ascending ModelID order. This gives disjoint donor identifiers between source
and target. Only lineages with N >= 10 contribute to fitting or testing.

For raw Chronos scores Y, the scientific effect is theta = mean(Y)/sd(Y).
Under the same independent Gaussian working model, test H0: theta >= 0
against H1: theta < 0. Computation uses X = -Y and delta = -theta > 0;
the figures retain the scientific Chronos sign.

For each held-out lineage, fit priors using only the other eligible lineages.
Gaussian, Cauchy, Laplace and Student-t5 use full-family location-scale MLE on
the standardized source estimates, one observation per source lineage.
Horseshoe uses the supplied constrained marginal MLE with noncentral-t
sampling distributions and source sample sizes, with scale >= 1e-4. The
fitted full families are then conditioned on delta > 0. Numerical multistart
optimization does not certify a global optimum.

Rouder is fixed Cauchy(0, 1/sqrt(2)) conditioned on delta > 0; it is not fitted.
The Lebesgue baseline is the actual improper flat measure on delta > 0,
with its separate extended-Ville calibration. There is no Ville result for
the improper baseline.

All settings are preserved: alpha = 0.05, monitoring from n0 = 2, 100 random
orders per lineage, base seed 20261009, quadrature order 48. The existing
seed rule includes gene, sex and lineage. Within each lineage, all procedures
use identical permutations. The same mixture evidence is used by eVille and
Ville for a given fitted prior; only calibration differs.

The table and raw-time plot require all evaluated random orders to reject;
they do not replace nonrejections with the final sample size. The full output
also retains the original general-purpose censored-summary columns.

## Interpretation

The RPS23 example was selected retrospectively to illustrate source-target
conflict and a difference in stopping times; see `SELECTION.md`. It is an
exploratory demonstration, not a new confirmatory biological finding or a
multiple-testing analysis of all screened genes. Calibration statements rely
on the working assumptions; donor deduplication alone does not establish
Gaussian sampling or independence.

Permutation means describe sensitivity to order conditional on the observed
data. They are not repeated-sampling power estimates. The different lineages
are not a controlled path of increasing prior conflict, so these figures do
not establish a dip and eventual return to the Lebesgue baseline.

Prostate is an outlier in its standardized mean, largely because its sample
standard deviation is larger. Its raw mean Chronos score is similar to those
of the source lineages. A prior's location is already on the standardized
effect scale; the displayed location-scale discrepancy is descriptive, not a
test statistic, and a Cauchy or horseshoe scale is not a standard deviation.

## Code map and reproducibility

- `cohort.py`, `effect_intervals.py`: filtering and fixed-sample effect intervals.
- `fitting.py`: the four MLE fits, horseshoe marginal MLE and fixed benchmarks.
- `engine.py`, `improper_flat.py`: evidence and calibration calculations.
- `sequential.py`: paired permutations, crossing boundaries and summaries.
- `run_analysis.py`: all held-out lineages and consolidated outputs.
- `validate.py`: numerical reference checks and independent output checks.
- `figures.py`: the same figure functions and style as the supplied package.

Only runner defaults/labels and the two Myeloid-specific presentation calls
were adapted. `method_preservation.json` records matching module hashes.
`environment.json` records versions used for the included run. See
`data/README.md` for the original DepMap release and file provenance.

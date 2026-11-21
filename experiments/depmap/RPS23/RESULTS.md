# Completed RPS23 results

Recorded sex: Male. Eligible cohort: **606 cell lines across 18 lineages**,
each with at least ten observations. Prostate has **10** observations;
its sources contain **596** observations across the other **17** lineages.
One ModelID per PatientID is retained before splitting.

## Prostate as the illustrative held-out target

Target standardized effect: -3.790417.
Full-family source Gaussian location: -8.294648; scale: 1.849646.
Descriptive discrepancy: 2.435185 Gaussian scale units.
This concerns the standardized effect, not a comparably large difference
in raw mean Chronos scores. Prostate has greater sample variability.

| Prior | Mean eVille | Mean Ville | Mean observations saved | Original eVille / Ville |
|---|---:|---:|---:|---|
| Gaussian | 3.51 | 3.94 | 0.43 | 2 / 6 |
| Cauchy | 3.58 | 3.91 | 0.33 | 6 / 6 |
| Laplace | 3.44 | 3.88 | 0.44 | 2 / 6 |
| Student-t5 | 3.47 | 3.92 | 0.45 | 2 / 6 |
| Horseshoe | 3.62 | 4.04 | 0.42 | 2 / 7 |
| Rouder half-Cauchy | 4.26 | 4.26 | 0.00 | 5 / 5 |
| Lebesgue | 3.56 | — | — | 4 / — |

The Gaussian eVille threshold is 14.152645, compared with 20 for Ville.
The Prostate horseshoe scale is 0.434824, above the fixed lower bound of 0.0001.
The same source-only marginal-likelihood fit as the uploaded package is used.

## Equal-lineage mean stopping times

These first average over 100 random orders within each lineage, then give
equal weight to the 18 eligible lineages. The illustrative ModelID order
is excluded.

| Fitted prior | Fitted eVille | Fitted Ville | Lebesgue eVille | Rouder eVille |
|---|---:|---:|---:|---:|
| Gaussian | 2.41 | 2.59 | 2.50 | 3.51 |
| Cauchy | 2.42 | 2.61 | 2.50 | 3.51 |
| Laplace | 2.41 | 2.59 | 2.50 | 3.51 |
| Student-t5 | 2.41 | 2.59 | 2.50 | 3.51 |
| Horseshoe | 2.42 | 2.60 | 2.50 | 3.51 |

Every procedure rejects in all 100 evaluated random orders in every
eligible lineage. Accordingly these are uncensored mean first rejection
times, not means capped at the available sample size. There are no Ville
values for the improper measure.

## What these results support

For the selected Prostate target, eVille has a smaller mean first rejection
time than Ville for every fitted prior. The raw stopping-time figure
also shows the other eligible targets and both fixed benchmarks.
Several fitted-prior eVille means are slightly below the Lebesgue mean,
whereas Cauchy and horseshoe are slightly above it. Thus this example
does not show a universal loss of the fitted prior advantage or a
controlled dip-and-return curve as conflict increases.

The fixed-order 2-versus-6 Gaussian crossing is an illustration; the
100-order means 3.51 versus 3.94 are the relevant order-averaged comparison.
These are conditional-on-data order summaries. The gene and illustrative
lineage were selected retrospectively; see SELECTION.md.

## Verification

The supplied validation routine passed. It checks likelihood ratios
against noncentral-t densities; horseshoe tail normalization; improper
mixture/calibration identities; source/target donor separation; cohort
eligibility; saved stopping decisions and summaries; and the example
evidence under quadrature refinement. Direct mixture values also agree
with boundary-based decisions for the fixed order and the first two
random orders in every lineage. Exact diagnostics are included in
`results/validation.json` and `results/boundary_checks.json`.

These checks verify the tested numerical cases and output consistency;
they do not prove all numerical errors are impossible or validate the
Gaussian and independence assumptions for the biological data.

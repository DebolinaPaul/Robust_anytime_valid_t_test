# Interpretation of the RPL38 experiment

This example shows modest improvements in stopping time from extended Ville
calibration. It does not establish a repeated-sample power curve or convergence
to the improper-flat baseline as prior conflict increases.

The numerical values below were reproduced with the final code using seed
20261009 and 100 random orders per target.

## Design

The supplied DepMap Public 24Q4 data retain 451 distinct PatientIDs after the
recorded-Female filter and the documented quality and deduplication rules.
There are 26 lineages: 17 main lineages with at least 10 observations, seven
smaller lineages with 2--8 observations, and two untestable singleton lineages.
All 24 testable lineages are held out in turn. Sources must contain at least
10 observations, and the complete target lineage is excluded from fitting.

The scientific test concerns negative mean Chronos score. Computation uses
the sign-reversed response and monitors from n=2 at alpha=0.05. Each target
uses 100 paired random observation orders, plus ascending ModelID order as
an illustration. The random-order summaries exclude this illustrative order.

Gaussian, Cauchy, Laplace and Student-t5 priors use source-only full-family
maximum likelihood fits, subsequently conditioned on the alternative
half-line. The exact horseshoe retains a documented median/IQR exception:
its unrestricted location-scale likelihood has no finite maximizer. There
is no separate IQR analysis. Rouder and improper-flat comparisons are fixed;
the improper prior has eVille calibration only.

## Main results

For the 17 main lineages, every evaluated random order rejects by the last
observation for every evaluated method. The useful comparison is therefore
when rejection occurs, rather than the final rejection fraction.

Across the five fitted priors and main lineages, mean paired savings range
from 0.06 to 0.30 observations. Their equal-lineage, equal-prior average is
0.169 observations. Most random orders produce a tie. With the Gaussian
prior, 63.35% of main-lineage random orders reject under eVille at n=2 and
93.53% reject by n=3. This proximity to the monitoring start leaves little
room for larger savings.

For the held-out Myeloid lineage, the original ModelID order gives eVille
rejection at n=3 and Ville rejection at n=7 for every fitted prior. This
four-observation saving is illustrative, not the typical saving over orders.
The random-order means are:

| Prior | Mean eVille stopping time | Mean Ville stopping time | Mean observations saved |
|---|---:|---:|---:|
| Gaussian | 3.61 | 3.82 | 0.21 |
| Cauchy | 3.48 | 3.76 | 0.28 |
| Laplace | 3.52 | 3.71 | 0.19 |
| Student-t5 | 3.51 | 3.78 | 0.27 |
| Horseshoe (quantile fit) | 3.54 | 3.71 | 0.17 |
| Rouder half-Cauchy | 4.30 | 4.30 | 0.00 |
| Improper flat | 3.36 | | |

## What can be said about conflict?

The descriptive conflict coordinate is

    abs(target standardized sample mean - source-fitted location) / fitted scale.

It uses an estimated target effect, not the unknown true mean. The scale is
the family's scale parameter, not necessarily a standard deviation; these
coordinates should be interpreted within each prior family. The coordinate
does not alter prior fitting, calibration, observation ordering or eligibility.

Myeloid has an estimated standardized effect of -3.65, while its fitted
Gaussian source prior has location -7.74 and scale 1.34: a discrepancy of
about 3.05 scale units. Its eVille procedure is faster than Ville, but the
flat baseline is still faster on average (3.36 versus 3.61 observations).
For each of the five fitted priors, eVille is slightly faster than flat in
the other 16 main lineages and slower than flat in Myeloid. This is consistent
with loss of the informative-prior benefit in this outlying target.

It is not evidence of a monotone general relationship: across the 17 main
lineages, the Spearman correlation between absolute Gaussian discrepancy
and mean paired saving is only about 0.10. Lineages also differ in effect
strength, sample size and observations. Conflict is not varied while these
other factors are held fixed.

To demonstrate a power-versus-conflict curve, use a controlled simulation
with repeated independent datasets and a specified prior perturbation.
A perturbation study on these fixed observations would instead describe
conditional sensitivity to prior choice and observation order. Convergence
to a flat baseline is specific to the prior family and limiting path; it
does not follow merely from an increasing discrepancy score.

## Small lineages and nonrejections

For nonrejections, restricted observation use is min(stopping time, N).
The savings bars separately distinguish equal rejections, earlier joint
rejections, eVille-only rejection, and neither method rejecting. Thus a
zero restricted saving is not presented as necessarily equal test outcomes.
For example, with N=2, Liver has no fitted-prior rejections, Pleura has both
procedures reject, and Vulva/Vagina has fitted eVille reject while Ville
does not; all three have zero restricted savings. Treat these groups as
small-sample sensitivity results, rather than pooling them into the main
17-lineage conclusion.

## Scope of the conclusions

Permutation percentages condition on the observed sample. They are not
repeated-sample power estimates or probabilities for new biological samples.
The normal sampling model and biological independence are working assumptions;
deduplicating PatientIDs does not prove independence. The gene, sex restriction
and illustrative target were chosen retrospectively. The analysis does not
provide selective-inference, familywise-error or false-discovery-rate control
for that search or the collection of lineage tests. Numerical checks support
the implementation, but do not certify every quadrature error or these
statistical assumptions.

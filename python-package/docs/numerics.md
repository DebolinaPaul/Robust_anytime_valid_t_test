# Numerical methods and calibration

## Statistical model and evidence

The model is an independent Gaussian sequence with unknown common variance.
The standardized effect is `(population_mean - mu) / sigma`. The package
uses scale-invariant likelihood ratios and mixes them against a fixed
effect-size prior. The resulting proper mixture is a nonnegative null
martingale in the invariant filtration, starting at observation 2.

The same mixture evidence is used for Ville and eVille. At level `alpha`,
Ville uses threshold `1/alpha`. The population extended Ville calibration
is `g(c) = E0[min(1, M2/c)]`. Its calibrated running value applies `g` to
the largest evidence observed so far. The boundary null has standardized
effect zero; for one-sided priors it is least favourable within the
corresponding composite null under the Gaussian model.

Full-line proper-prior evidence uses the simulation notebook's normal-scale
mixture calculations, including its Gaussian component formula, stable
logarithmic computations and quadrature refinement. One-sided evidence
integrates over the selected effect half-line and normalizes a proper
prior by its half-line probability. A one-sided calibration must use this
one-sided evidence; it is not the two-sided calibration reused unchanged.

The horseshoe is the genuine normal--half-Cauchy mixture. Its density has
an integrable infinite spike at its location. Tail normalization and
integration near that spike need separate numerical treatment.

## Three different errors

1. **Evidence integration.** Refinement and tail checks assess numerical
   integration of the prior mixture. These are floating-point calculations,
   not interval-arithmetic certificates.
2. **Evaluation of the empirical calibration.** For the fixed sample of
   `B` null statistics, adaptive convex bounds bracket the finite average
   of `min(1, M2_null/c)`, assuming accurate evidence evaluations. The
   `mc_lower` and `mc_upper` columns describe this numerical bracket.
3. **Sampling the null distribution.** The finite average differs from
   the population expectation. `mcse` describes pointwise Monte Carlo
   uncertainty. It is not a simultaneous confidence band and does not
   restore an exact anytime-valid guarantee.

Reducing numerical tolerances addresses the first two errors, not the
third. Numerical bracket endpoints must not be presented as confidence
limits for the population calibration.

The default empirical upper value is capped by `min(1, 1/c)`, as in the
notebook. This enforces the comparison with ordinary Ville for the same
proper mixture; it does not compensate for possible underestimation of
the population extended Ville calibration by a finite Monte Carlo sample.

The deterministic option evaluates the population expectation through
analytic formulas or numerical integration. It avoids Monte Carlo sampling
error, while retaining the usual numerical limitations of floating-point
evidence and integration calculations. It is not inherently less
statistically defensible than Monte Carlo sampling.

## Reproducibility

The bundled sample contains 100,000 t statistics from independent pairs of
standard-normal observations, generated with NumPy seed `20260510` in the
paper's simulation notebook. Bundling the actual statistics makes agreement
independent of differences between R and Python random-number generators.

The standard settings are absolute calibration tolerance `1e-5`, trajectory
relative tolerance `1e-3`, and the notebook's proper-prior parameterizations.
Matching an existing calculation requires matching all of these settings
as well as the target observations, their order, and the prior.

The bundled sample can be reused across priors or target datasets; it must
remain independent of the target observations. This reuse gives comparable
Monte Carlo calculations but does not create independent calibration
estimates for each test.

## Lebesgue benchmarks

Lebesgue measure is an improper mixing measure, not a probability density.
Its mixture evidence is finite almost surely from observation 2 under the
Gaussian model but has infinite null expectation. The package applies its
separate analytic extended Ville calibration, using the full-line or
half-line formula as appropriate. Ordinary reciprocal Ville calibration
is unavailable for this measure.

Multiplying the mixing measure by a positive constant also rescales its
evidence and calibration, leaving the calibrated result unchanged. The
benchmark need not be the most powerful choice at every sample size or
significance level.

## Interpreting exported trajectories

`log_evidence` describes the evidence at the current prefix;
`running_log_evidence` retains its running maximum. `p_eville` and
`p_ville` are running calibrated values. Under the default method,
`p_eville` is the Monte Carlo approximation described above.

`calibration_lower` and `calibration_gap` describe the numerical bounds
after the Ville cap and running-minimum operation. The `mc_*` fields retain
the raw finite-sample calibration diagnostics at the current running
evidence threshold, before that cap. These fields therefore have distinct
interpretations and need not coincide.

The first rejection time is the earliest supplied prefix with a running
value at most `alpha`. No rejection within the observed sequence is not
a statement about all future observations. Changing the prior, direction,
or observation order after inspecting the target is not covered by the
fixed-prior guarantee.

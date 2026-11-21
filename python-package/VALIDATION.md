# Validation

The package was checked on 21 November 2026 using Python 3.12.14,
NumPy 2.3.5, SciPy 1.17.0, pandas 2.2.3 and Matplotlib 3.10.8 on Linux.
The source distribution and wheel build successfully. The installed wheel
loads its bundled null sample outside the source directory.

The test suite has 49 passing tests, covering:

- Densities for all five proper prior families.
- Mixture evidence against independently generated notebook references.
- One-sided integration against SciPy noncentral-t densities.
- Half-line recombination, reflection and positive-rescaling invariance.
- Extreme horseshoe tail normalization and an integration regression near
  the horseshoe density spike.
- The shared 100,000-draw Monte Carlo calibration, and numerical brackets
  checked against direct finite-sample averages.
- Analytic Gaussian and full/half-Lebesgue calibration.
- General deterministic calibration against independent reference values.
- First rejection times, summary tables, plot generation and invalid inputs.

Additional numerical checks compared 30 deterministic calibrations across
all five families and all three alternatives with independent adaptive
integration of the null-angle expectation. Those integrals lie within the
reported numerical bounds, up to the reference integration tolerance.

The two-sided evidence engine is the numerical engine from
`experiments/simulations/Code_simulation.ipynb`. The Monte Carlo engine
uses the same empirical calibration rule. The notebook SHA256 is
`799a35a2392812cacd12491bc40b3c9c1de68f2039ffa550d5b83d764e7dfccc`.
The bundled sample and its provenance are included in `src/eville/data/`.
Tests use reference CSVs generated separately from this package.

Both examples run successfully and export CSV and PDF files. For the
100-observation stream used in `examples/compare_priors.py`, the default
settings give first rejection times 48, 47, 47, 47 and 48 for Gaussian,
Cauchy, Student-t5, Laplace and horseshoe, respectively, for both methods.
This example demonstrates the interface and agreement with the reference
calculations; it is not an example of a strict stopping-time improvement.

A separate precision check reproduced all 15 Monte Carlo reference rows
when extended-precision accumulations were restricted to ordinary float64,
as on some Mac builds. The maximum difference was below `7e-16`.
This does not replace running the checks on your Mac:

```bash
python -m pip install -e ".[dev]"
python -m pytest
python examples/quick_start.py
python -m build
```

Tests and numerical comparisons do not prove mathematical validity. The
theory assumes the stated Gaussian model and a prior fixed independently
of the target observations. The default finite Monte Carlo calibration is
an approximation to the population calibration; its numerical bracket is
not a bound on Monte Carlo sampling error. Floating-point quadrature
checks are not interval-arithmetic certificates.

"""Regenerate R reference fixtures from the repository's simulation notebook.

Run from a checkout with the notebook's Python dependencies installed:
    python r-package/tools/export_all_priors_reference.py

Only definition cells are executed; simulations and figure cells are not run.
The bundled null sample preserves the notebook's actual NumPy draws, rather
than assuming that equal seed integers generate equal R and NumPy streams.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import platform
import sys
import types

import numpy as np
import scipy
from scipy import stats


def write_csv(path, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    package = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--notebook", type=Path,
                        default=package.parent / "experiments/simulations/Code_simulation.ipynb")
    args = parser.parse_args()
    notebook_bytes = args.notebook.read_bytes()
    notebook = json.loads(notebook_bytes)
    module = types.ModuleType("eville_reference_export")
    sys.modules[module.__name__] = module
    namespace = module.__dict__
    for index in (4, 6, 8):
        cell = notebook["cells"][index]
        if cell["cell_type"] != "code":
            raise ValueError("The notebook cell layout changed; review definition-cell indices.")
        exec(compile("".join(cell["source"]), f"{args.notebook.name}:cell{index}", "exec"), namespace)

    config = namespace["simulation_config"]("paper")
    seed = config["seed_calibration"]
    B = config["calibration_replicates"]
    if seed != 20260510 or B != 100000:
        raise ValueError("Notebook calibration settings changed; update the bundled sample's name.")
    draws = np.random.default_rng(seed).normal(0.0, 1.0, (B, 2))
    t_null = namespace["t_statistic"](draws)
    extdata = package / "inst/extdata"
    extdata.mkdir(parents=True, exist_ok=True)
    null_path = extdata / "python_null_t_seed20260510.csv.gz"
    text = "t_null\n" + "".join(f"{value:.17g}\n" for value in t_null)
    with null_path.open("wb") as stream:
        with gzip.GzipFile(filename="", mode="wb", fileobj=stream, mtime=0) as compressed:
            compressed.write(text.encode("ascii"))

    fixtures = package / "tests/testthat/fixtures"
    fixtures.mkdir(parents=True, exist_ok=True)
    families = {"gaussian": "Gaussian", "cauchy": "Cauchy",
                "student_t5": "Student-t5", "laplace": "Laplace",
                "horseshoe": "Horseshoe"}
    cases = [
        (2, -3.2, -8, 0.7), (2, 1.2, 0, 1), (2, 15, 12, 1),
        (5, -2.4, 3, 0.5), (5, 0.7, -12, 1),
        (20, 2.8, 0.5, 1), (20, -4, -12, 2),
        (50, 3.5, 0, 1), (50, -4.5, 12, 1),
        (50, 5.1, 50, 0.7), (2, -0.8, 50, 1),
    ]
    evidence = []
    for family, python_family in families.items():
        for n, t, location, scale in cases:
            prior = namespace["Prior"](python_family, location, scale)
            diagnostic = namespace["log_mixture_with_diagnostics"](t, n, prior)
            evidence.append(dict(family=family, n=n, t=t, location=location,
                                 scale=scale, log_evidence=diagnostic.log_value,
                                 relative_refinement_change=diagnostic.relative_refinement_change,
                                 relative_omitted_tail_bound=diagnostic.relative_omitted_tail_bound))
    write_csv(fixtures / "all_priors_evidence.csv", evidence)

    densities = []
    for family, python_family in families.items():
        location, scale = -2.0, 1.3
        for z in (-100.0, -8.0, -0.01, 0.01, 1.0, 8.0, 100.0):
            theta = location + scale * z
            if family == "gaussian":
                lp = stats.norm.logpdf(theta, loc=location, scale=scale)
            elif family == "cauchy":
                lp = stats.cauchy.logpdf(theta, loc=location, scale=scale)
            elif family == "student_t5":
                lp = stats.t.logpdf(theta, df=5, loc=location, scale=scale)
            elif family == "laplace":
                lp = stats.laplace.logpdf(theta, loc=location, scale=scale)
            else:
                lp = namespace["horseshoe_logpdf"](theta, location, scale)
            densities.append(dict(family=family, theta=theta, location=location,
                                  scale=scale, log_density=float(lp)))
    write_csv(fixtures / "all_priors_densities.csv", densities)

    # Direct evaluation at every sampled t verifies the empirical mean without
    # relying on the adaptive affine approximation used by either language.
    small_B = 257
    direct = []
    for family, python_family in families.items():
        scale = 1 / np.sqrt(2) if family == "cauchy" else 1.0
        location = 3.0
        prior = namespace["Prior"](python_family, location, scale)
        logs = np.array([namespace["log_mixture"](float(t), 2, prior)
                         for t in t_null[:small_B]])
        for log_c in (-3.0, 0.0, np.log(20), 7.0):
            values = np.exp(np.minimum(0.0, logs - log_c))
            estimate = float(values.mean())
            mcse = float(values.std(ddof=1) / np.sqrt(small_B))
            ville = float(np.exp(min(0, -log_c)))
            direct.append(dict(family=family, location=location, scale=scale,
                               B=small_B, log_c=log_c, empirical_mean=estimate,
                               direct_mcse=mcse, ville_bound=ville,
                               used_if_exact=min(estimate, ville)))
    write_csv(fixtures / "all_priors_mc_direct.csv", direct)

    # Paper-sized, identical-null-sample reference for the adaptive algorithm.
    null_sample = namespace["NullSample"](t_null)
    adaptive = []
    tolerance = config["calibration_numerical_tolerance"]
    for family, python_family in families.items():
        scale = 1 / np.sqrt(2) if family == "cauchy" else 1.0
        location = 2.0 if family == "gaussian" else 8.0
        prior = namespace["Prior"](python_family, location, scale)
        cal = namespace["MCCalibration"](prior, null_sample, tolerance=tolerance,
                                          max_intervals=config["calibration_max_intervals"])
        for log_c in (0.0, np.log(20), 7.0):
            result = cal.evaluate(log_c)
            adaptive.append(dict(family=family, location=location, scale=scale,
                                 B=B, log_c=log_c, tolerance=tolerance,
                                 estimate=result.estimate, lower=result.lower,
                                 upper=result.upper, gap=result.gap,
                                 mcse=result.mcse, used=result.used,
                                 ville_bound=result.ville_bound,
                                 intervals=result.intervals))
    write_csv(fixtures / "all_priors_mc_100000.csv", adaptive)

    metadata = dict(
        source="experiments/simulations/Code_simulation.ipynb",
        source_sha256=hashlib.sha256(notebook_bytes).hexdigest(),
        definition_cells=[4, 6, 8],
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
        calibration_seed=seed, calibration_replicates=B,
        null_generation="numpy.random.default_rng(seed).normal(0,1,(B,2)); notebook t_statistic",
        null_order="Unsorted, original NumPy draw order; calibrators sort a copy.",
        null_csv_sha256=hashlib.sha256(text.encode("ascii")).hexdigest(),
        null_gzip_sha256=hashlib.sha256(null_path.read_bytes()).hexdigest(),
        first_five_null_t=t_null[:5].tolist(),
        paper_calibration_tolerance=tolerance,
        evidence_rows=len(evidence), density_rows=len(densities),
        direct_calibration_rows=len(direct),
        adaptive_calibration_rows=len(adaptive),
        note="Empirical numerical brackets do not bound Monte Carlo error or establish exact population anytime validity."
    )
    (fixtures / "all_priors_reference_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

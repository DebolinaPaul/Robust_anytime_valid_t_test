#!/usr/bin/env python3
"""Export TWO-SIDED regression fixtures from the existing RPL38 Python engine.

Run from the eVille R-package root, for example:
  python tools/export_python_reference.py \
    --engine-dir /path/to/RPL38_final

Requires numpy and scipy, and engine.py plus improper_flat.py in --engine-dir.
No target/source data files are used. The input vector is deliberately fixed.

The RPL38 engine integrates over POSITIVE effects. For a centered symmetric
proper prior, its full two-sided mixture is the average of the positive-half
mixtures at t and -t. For Lebesgue measure it is the SUM, with no division by 2.
The RPL38 one-sided calibration is deliberately NOT used as a reference for
two-sided calibration.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import math
import platform
from pathlib import Path
import sys

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import logsumexp
from scipy.stats import nct, t as student_t


X = np.array([-0.4, 0.9, 0.2, 1.1, -0.7, 0.5,
              0.3, 1.4, -0.2, 0.8, 0.1, 0.6])


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def full_log_mixture(engine, flat, t, n, family, scale, order):
    if family == "lebesgue":
        return float(np.logaddexp(flat.log_flat_mixture(t, n),
                                 flat.log_flat_mixture(-t, n)))
    name = "Gaussian" if family == "gaussian" else "Cauchy"
    prior = engine.Prior(name, 0.0, scale)
    halves = [engine.log_mixture(sign * t, n, prior,
                                quadrature_order=order)
              for sign in (1, -1)]
    return float(logsumexp(halves) - math.log(2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path,
                        default=Path("tests/testthat/fixtures"))
    args = parser.parse_args()
    engine_dir = args.engine_dir.expanduser().resolve(strict=True)
    paths = {name: engine_dir / name for name in ("engine.py", "improper_flat.py")}
    for path in paths.values():
        if not path.is_file():
            parser.error(f"Missing {path}")
    sys.path.insert(0, str(engine_dir))
    engine = importlib.import_module("engine")
    flat = importlib.import_module("improper_flat")
    for module, name in ((engine, "engine.py"), (flat, "improper_flat.py")):
        if Path(module.__file__).resolve() != paths[name]:
            raise RuntimeError(f"Imported unexpected {name}")

    out = args.output_dir.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "python_observations.csv",
              [{"i": i + 1, "x": float(value)} for i, value in enumerate(X)])

    rows, refinements, analytic_errors = [], [], []
    settings = [("gaussian", 1.0), ("gaussian", 3.0),
                ("cauchy", 1 / math.sqrt(2)), ("lebesgue", None)]
    cases = [("prefix", n, math.sqrt(n) * X[:n].mean() / X[:n].std(ddof=1))
             for n in range(2, len(X) + 1)]
    cases += [("zero_t", n, 0.0) for n in (2, 5, 20)]
    for case, n, t in cases:
        d = n * (n - 1) / ((n - 1) + t * t)
        for family, scale in settings:
            low = full_log_mixture(engine, flat, t, n, family, scale, 32)
            high = full_log_mixture(engine, flat, t, n, family, scale, 64)
            refinements.append(abs(high - low))
            if family == "gaussian":
                exact = ((n - 1) / 2 * math.log1p(n * scale * scale)
                         - n / 2 * math.log1p(d * scale * scale))
                analytic_errors.append(abs(high - exact))
            elif family == "lebesgue":
                exact = 0.5 * math.log(2 * math.pi / n) + n / 2 * math.log(n / d)
                analytic_errors.append(abs(high - exact))
            else:
                # Independent Cauchy normal-scale-mixture integral. This
                # oracle integrates an unbounded interval, not the old
                # engine's finite positive-effect grid.
                def integrand(w):
                    return math.exp(-w / 2 + (n - 1) / 2 * math.log(w + n * scale**2)
                                    - n / 2 * math.log(w + d * scale**2)
                                    - 0.5 * math.log(2 * math.pi))
                value, error = quad(integrand, 0, np.inf,
                                    epsabs=2e-12, epsrel=2e-12, limit=300)
                if error > 1e-10 * max(1.0, value):
                    raise ArithmeticError("Independent Cauchy check failed quadrature tolerance")
                analytic_errors.append(abs(high - math.log(value)))
            rows.append({"case": case, "n": n, "t": float(t),
                         "prior": family, "scale": "" if scale is None else scale,
                         "log_m": high})

    # Moderate values only: a separate likelihood-ratio check using scipy's
    # noncentral-t implementation, with the correct sqrt(n)*theta convention.
    lr_rows, lr_errors = [], []
    for n in (2, 5, 20):
        for t in (-2.0, 0.0, 2.0):
            for theta in (0.0, 0.3, 1.0):
                actual = float(engine.log_likelihood_ratio(t, n, theta))
                oracle = float(nct.logpdf(t, n - 1, math.sqrt(n) * theta)
                               - student_t.logpdf(t, n - 1))
                lr_errors.append(abs(actual - oracle))
                lr_rows.append({"n": n, "t": t, "theta": theta,
                                "log_lr": actual})

    report = {
        "scope": "two-sided centered priors; evidence only; n0 not involved",
        "source_files_sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest()
                                for name, path in paths.items()},
        "versions": {"python": platform.python_version(),
                     "numpy": np.__version__, "scipy": scipy.__version__},
        "quadrature_orders": [32, 64],
        "max_log_m_order_difference": max(refinements),
        "max_log_m_independent_formula_difference": max(analytic_errors),
        "max_log_lr_nct_difference": max(lr_errors),
        "mixture_rows": len(rows),
        "likelihood_ratio_rows": len(lr_rows),
        "proper_symmetrization": "logaddexp(log_Mplus(t),log_Mplus(-t))-log(2)",
        "lebesgue_symmetrization": "logaddexp(log_Mplus(t),log_Mplus(-t))",
    }
    if max(refinements + analytic_errors + lr_errors) > 1e-8:
        raise ArithmeticError(json.dumps(report, indent=2))
    write_csv(out / "python_evidence.csv", rows)
    write_csv(out / "python_likelihood_ratio.csv", lr_rows)
    (out / "python_reference_metadata.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

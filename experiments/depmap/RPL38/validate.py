#!/usr/bin/env python3
"""Independent numerical and saved-output checks for the final RPL38 analysis.

These checks establish numerical agreement and bookkeeping consistency for
the tested cases. They do not certify quadrature error, Gaussian sampling,
biological independence, or validity after retrospective example selection.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.stats import nct, norm, t as student_t

import engine
from improper_flat import FlatPrior, flat_calibration, flat_log_threshold, log_flat_mixture


FLAT = "Improper flat"


def require(condition, message):
    if not bool(condition):
        raise AssertionError(message)


def close(actual, expected, context, rtol=1e-9, atol=1e-10):
    if expected is None or pd.isna(expected):
        require(pd.isna(actual), f"{context}: expected a blank value")
    elif isinstance(expected, str):
        require(actual == expected, f"{context}: {actual!r} != {expected!r}")
    else:
        require(np.isclose(float(actual), float(expected), rtol=rtol, atol=atol),
                f"{context}: {actual!r} != {expected!r}")


def validate_core():
    errors = []
    for n in (2, 5, 10, 50):
        for t in (-3.0, -0.5, 0.0, 0.5, 3.0):
            for delta in (0.0, 0.2, 0.6):
                reference = nct.logpdf(t, n - 1, math.sqrt(n) * delta) - student_t.logpdf(t, n - 1)
                actual = float(engine.log_likelihood_ratio(t, n, delta))
                errors.append(abs(actual - reference))
    require(max(errors) < 1e-8, "Likelihood ratio differs from the independent noncentral-t reference")

    tail_errors = []
    for z in (1.0, 100.0, 1e4, 1e8):
        # Independent mixture representation: integrate over the normal
        # variate, rather than the survival transform used in the engine.
        value, error = quad(
            lambda u: (2 / math.pi) * norm.pdf(u) * z * math.atan(u / z),
            0.0, np.inf, epsabs=1e-13, epsrel=1e-12, limit=200,
        )
        reference = value / z
        actual = engine._hs_sf_positive(z)
        tail_errors.append(abs(actual / reference - 1))
        require(error < 1e-10, "Independent horseshoe-tail integration failed")
    require(max(tail_errors) < 1e-9, "Horseshoe tail normalization failed")

    prior = engine.Prior("Horseshoe", -10000.0, 1.0)
    reference, _ = quad(
        lambda d: math.exp(-d * d + float(prior.logpdf(d))),
        0.0, np.inf, epsabs=1e-13, epsrel=1e-11, limit=200,
    )
    hs_log_error = abs(engine.log_mixture(0.0, 2, prior, 96) - math.log(reference))
    require(hs_log_error < 1e-7, "Extreme horseshoe mixture normalization failed")

    flat_errors = []
    for n in (2, 3, 10, 30):
        for t in (-5.0, -1.0, 0.0, 1.0, 5.0):
            expected = log_flat_mixture(t, n)
            center = max(0.0, t / math.sqrt(n))
            width = math.hypot(t, math.sqrt(n - 1)) / math.sqrt(n * (n - 1))
            points = sorted(set([0.0, center, center + width, math.inf]))
            integral = sum(quad(
                lambda d: math.exp(float(engine.log_likelihood_ratio(t, n, d)) - expected),
                a, b, epsabs=1e-10, epsrel=1e-10, limit=200,
            )[0] for a, b in zip(points[:-1], points[1:]))
            flat_errors.append(abs(integral - 1))
    require(max(flat_errors) < 2e-8, "Flat mixture fails independent likelihood integration")

    inverse_errors = []
    for alpha in (1e-300, 1e-100, 1e-16, 0.001, 0.05, 0.5, 0.99):
        p = flat_calibration(flat_log_threshold(alpha))
        inverse_errors.append(abs(p / alpha - 1))
    require(max(inverse_errors) < 1e-10, "Flat calibration inversion failed")

    near_zero_errors = []
    for n, t in ((10, 1e-8), (1000, 1e-7), (50, -1e-8)):
        reference = (0.5 * math.log(2 * math.pi / n)
                     + 0.5 * n * math.log1p(t * t / (n - 1))
                     + student_t.logcdf(t * math.sqrt(n / (n - 1)), n))
        near_zero_errors.append(abs(log_flat_mixture(t, n) - reference))
    require(max(near_zero_errors) < 1e-10, "Near-zero flat mixture lost precision")

    invalid_cases = [(0.0, 2.9, 0.1), (np.nan, 2, 0.1), (0.0, 2, np.nan), (0.0, 2, -0.1)]
    for args in invalid_cases:
        try:
            engine.log_likelihood_ratio(*args)
        except (ValueError, FloatingPointError):
            pass
        else:
            raise AssertionError(f"Invalid likelihood input accepted: {args}")
    return {
        "likelihood_reference_points": len(errors),
        "maximum_likelihood_log_error": max(errors),
        "maximum_horseshoe_survival_relative_error": max(tail_errors),
        "extreme_horseshoe_log_mixture_error": hs_log_error,
        "flat_mixture_reference_points": len(flat_errors),
        "maximum_flat_mixture_relative_error": max(flat_errors),
        "maximum_flat_calibration_inverse_relative_error": max(inverse_errors),
        "maximum_near_zero_flat_log_error": max(near_zero_errors),
        "invalid_likelihood_inputs_rejected": len(invalid_cases),
    }


def independent_summary(tau, N):
    tau = np.asarray(tau, dtype=float)
    observed = np.isfinite(tau)
    finite = tau[observed]
    require(np.all((finite >= 2) & (finite <= N) & (finite == np.floor(finite))),
            "Invalid finite first-rejection time")
    capped = np.where(observed, tau, N)
    median = np.sort(np.where(observed, tau, np.inf))[math.ceil(len(tau) / 2) - 1]
    return {
        "random_permutations": len(tau),
        "rejection_count": int(observed.sum()),
        "rejection_fraction": float(observed.mean()),
        "mean_capped_at_N": float(capped.mean()),
        "median_capped_at_N": float(np.median(capped)),
        "mean_capped_at_N_over_N": float(capped.mean() / N),
        "median_capped_at_N_over_N": float(np.median(capped) / N),
        "unconditional_median_rejection_time": int(median) if np.isfinite(median) else None,
        "unconditional_median_status": "observed" if np.isfinite(median) else "NR",
        "mean_first_rejection_time_if_all_orders_reject": float(tau.mean()) if observed.all() else None,
        "conditional_on_rejection_mean": float(finite.mean()) if len(finite) else None,
        "conditional_on_rejection_median": float(np.median(finite)) if len(finite) else None,
    }


def independent_pairs(e, v, N):
    e = np.where(np.isfinite(e), e, np.inf)
    v = np.where(np.isfinite(v), v, np.inf)
    require(np.all(e <= v), "eVille rejects after Ville")
    both = np.isfinite(e) & np.isfinite(v)
    difference = v[both] - e[both]
    restricted = np.minimum(v, N) - np.minimum(e, N)
    return {
        "random_permutations": len(e),
        "eville_earlier_fraction": float(np.mean(e < v)),
        "both_reject_fraction": float(np.mean(both)),
        "eville_only_rejects_fraction": float(np.mean(np.isfinite(e) & ~np.isfinite(v))),
        "neither_rejects_fraction": float(np.mean(~np.isfinite(e) & ~np.isfinite(v))),
        "earlier_and_both_reject_fraction": float(np.mean((e < v) & both)),
        "mean_restricted_observations_saved": float(restricted.mean()),
        "median_restricted_observations_saved": float(np.median(restricted)),
        "mean_saved_conditional_on_both_rejecting": float(difference.mean()) if len(difference) else None,
        "median_saved_conditional_on_both_rejecting": float(np.median(difference)) if len(difference) else None,
    }


def validate_results(results):
    """Validate consolidated outputs without using their summary functions."""
    results = Path(results)
    settings = json.loads((results / "settings.json").read_text())
    observations = pd.read_csv(results / "observations.csv.gz")
    cohort = pd.read_csv(results / "cohort.csv")
    stopping = pd.read_csv(results / "stopping_times.csv.gz")
    summary = pd.read_csv(results / "summary.csv")
    fits = pd.read_csv(results / "prior_fits.csv")
    boundaries = pd.read_csv(results / "critical_boundaries.csv.gz")
    example = pd.read_csv(results / "example_paths.csv")
    B = int(settings["permutations"])
    alpha = float(settings["alpha"])
    priors = set(engine.FAMILIES) | {"Rouder half-Cauchy", FLAT}
    require(B >= 1 and settings["n0"] == 2, "Unsupported monitoring configuration")
    require(not observations.ModelID.duplicated().any(), "Duplicate retained ModelID")
    require(not observations.PatientID.duplicated().any(), "Duplicate retained PatientID")
    require(observations.PatientID.notna().all(), "Missing retained PatientID")
    require(np.isfinite(observations.Chronos).all(), "Nonfinite retained score")
    require(observations.Sex.eq(settings["sex"]).all(), "Recorded-sex filter differs")
    require(len(observations) == settings["observations_retained"], "Retained sample count differs")
    require(set(cohort.target_lineage) == set(observations.OncotreeLineage), "Cohort inventory incomplete")
    require(not cohort.target_lineage.duplicated().any(), "Duplicated cohort lineage")
    require(not summary.duplicated(["target_lineage", "prior"]).any(), "Duplicated summary key")
    require(not fits.duplicated(["target_lineage", "prior"]).any(), "Duplicated prior-fit key")
    require(not stopping.duplicated(["target_lineage", "prior", "permutation"]).any(),
            "Duplicated stopping-time key")

    counts = observations.groupby("OncotreeLineage").size()
    sds = observations.groupby("OncotreeLineage").Chronos.std(ddof=1)
    eligible = {name for name in counts.index
                if counts[name] >= settings["min_target_n"] and np.isfinite(sds[name]) and sds[name] > 0}
    sources_allowed = {name for name in counts.index
                       if counts[name] >= settings["min_source_n"] and np.isfinite(sds[name]) and sds[name] > 0}
    require(set(settings["tested_lineages"]) == eligible, "Tested lineages differ from eligibility rule")
    require(set(summary.target_lineage) == eligible, "Missing or extra tested lineage in summaries")
    require(set(stopping.target_lineage) == eligible, "Missing or extra tested lineage in raw records")
    require(set(fits.target_lineage) == eligible, "Missing or extra tested lineage in fits")
    main = {name for name in eligible if counts[name] >= settings["main_min_n"]}
    require(set(settings["main_lineages"]) == main, "Main cohort classification differs")
    require(set(settings["small_lineages"]) == eligible - main, "Small cohort classification differs")

    outcome_names = ["outcome_same_rejection_fraction", "outcome_save1_fraction",
                     "outcome_save2_fraction", "outcome_save3plus_fraction",
                     "eville_only_rejects_fraction", "neither_rejects_fraction"]
    paired_names = list(independent_pairs(np.array([2.]), np.array([2.]), 2))
    paired_names.remove("random_permutations")
    checked_decisions = checked_summaries = 0
    boundary_rows = 0
    example_ts = None
    for lineage in sorted(eligible):
        rows = observations.loc[observations.OncotreeLineage.eq(lineage)].sort_values("ModelID", kind="stable")
        x = -rows.Chronos.to_numpy(float)
        N = len(x)
        source_names = sources_allowed - {lineage}
        source_rows = observations.loc[observations.OncotreeLineage.isin(source_names)]
        require(not (set(rows.PatientID) & set(source_rows.PatientID)), f"Source/target donor overlap: {lineage}")
        require(lineage not in set(source_rows.OncotreeLineage), f"Target lineage leaked into source: {lineage}")

        # Reconstruct the published seed directly, without calling the runner.
        key = f"{int(settings['seed'])}\0{settings['gene']}\0{settings['sex']}\0{lineage}".encode()
        seed = int.from_bytes(hashlib.sha256(key).digest()[:8], "little")
        rng = np.random.default_rng(seed)
        indices = np.vstack([np.arange(N)] + [rng.permutation(N) for _ in range(B)])
        sample = x[indices]
        ns = np.arange(2, N + 1)
        # Direct centered standard deviations are independent of Welford updates.
        ts = np.stack([np.sqrt(float(n)) * sample[:, :n].mean(axis=1)
                       / sample[:, :n].std(axis=1, ddof=1) for n in ns], axis=1)
        require(np.isfinite(ts).all() and np.max(np.abs(ts)) < math.sinh(32),
                f"Invalid t statistics or boundary-search range: {lineage}")
        if lineage == settings["example_target"]:
            example_ts = ts[0]

        require(set(summary.loc[summary.target_lineage.eq(lineage), "prior"]) == priors,
                f"Incorrect prior set: {lineage}")
        flat_group = stopping.loc[stopping.target_lineage.eq(lineage) & stopping.prior.eq(FLAT)]
        flat_random = flat_group.loc[flat_group.permutation.gt(0), "eville_rejection"].to_numpy(float)
        flat_mean = float(np.where(np.isfinite(flat_random), flat_random, N).mean())
        for prior in sorted(priors):
            context = f"{lineage}/{prior}"
            raw = stopping.loc[stopping.target_lineage.eq(lineage) & stopping.prior.eq(prior)].sort_values("permutation")
            require(np.array_equal(raw.permutation.to_numpy(), np.arange(B + 1)), f"Incorrect paired orders: {context}")
            record = summary.loc[summary.target_lineage.eq(lineage) & summary.prior.eq(prior)].iloc[0]
            fit = fits.loc[fits.target_lineage.eq(lineage) & fits.prior.eq(prior)].iloc[0]
            close(record.target_n, N, context + "/target_n")
            close(record.source_n, len(source_rows), context + "/source_n")
            close(record.source_lineages, len(source_names), context + "/source_lineages")
            close(record.cohort, "main" if lineage in main else "small", context + "/cohort")
            close(record.theta_hat, x.mean() / -x.std(ddof=1), context + "/theta_hat")
            methods = ["eville"] if prior == FLAT else ["eville", "ville"]

            for prefix in methods:
                times = raw[prefix + "_rejection"].to_numpy(float)
                require(np.array_equal(raw[prefix + "_rejected"].eq(True).to_numpy(), np.isfinite(times)),
                        f"Rejection indicators disagree: {context}/{prefix}")
                expected = independent_summary(times[1:], N)
                for name, value in expected.items():
                    close(record[prefix + "_" + name], value, f"{context}/{prefix}_{name}")
                close(record["original_" + prefix + "_rejection"], times[0], context + "/original")
                checked_summaries += 1

                b = boundaries.loc[boundaries.target_lineage.eq(lineage)
                                   & boundaries.prior.eq(prior)
                                   & boundaries.method.eq("eVille" if prefix == "eville" else "Ville")].sort_values("n")
                require(np.array_equal(b.n.to_numpy(), ns), f"Boundary sample sizes differ: {context}/{prefix}")
                values = b.critical_t.to_numpy(float)
                expected_status = np.where(np.isnan(values), "not_needed_after_all_orders_rejected",
                    np.where(np.isfinite(values), "finite_root",
                             np.where(values > 0, "above_searched_range", "below_searched_range")))
                require(np.array_equal(b.boundary_status.to_numpy(), expected_status), f"Boundary status differs: {context}")
                crossed = ts >= values[None, :]
                direct = np.where(crossed.any(axis=1), crossed.argmax(axis=1) + 2, N + 1)
                expected_times = np.where(np.isfinite(times), times, N + 1)
                require(np.array_equal(direct, expected_times), f"Saved crossing differs from independent prefix statistics: {context}/{prefix}")
                omitted = np.flatnonzero(np.isnan(values))
                if len(omitted):
                    first = int(omitted[0])
                    require(np.isnan(values[first:]).all(), f"Non-trailing omitted boundary: {context}")
                    require(np.all(direct < ns[first]), f"Boundary omitted before all orders rejected: {context}")
                checked_decisions += len(times)
                boundary_rows += len(b)

            e = raw.eville_rejection.iloc[1:].to_numpy(float)
            close(record.flat_mean_capped_at_N, flat_mean, context + "/flat_mean")
            e_mean = independent_summary(e, N)["mean_capped_at_N"]
            close(record.eville_mean_minus_flat, e_mean - flat_mean, context + "/eville_minus_flat")

            if prior == FLAT:
                blank = [c for c in summary if c.startswith("ville_") or c == "original_ville_rejection"]
                blank += paired_names + outcome_names + ["conflict_signed_location_scale", "conflict_abs_location_scale"]
                require(record[list(dict.fromkeys(blank))].isna().all(), f"Improper-prior comparison fields must be blank: {context}")
                require(raw[["ville_rejection", "ville_rejected"]].isna().all().all(), f"Improper Ville raw fields must be blank: {context}")
                require(pd.isna(fit.ville_threshold), f"Improper Ville threshold must be blank: {context}")
                require(not ((boundaries.target_lineage.eq(lineage)) & boundaries.prior.eq(FLAT)
                             & boundaries.method.eq("Ville")).any(), f"Improper Ville boundary exists: {context}")
                close(flat_calibration(fit.log_eville_threshold), alpha, context + "/flat_calibration")
            else:
                v = raw.ville_rejection.iloc[1:].to_numpy(float)
                for name, value in independent_pairs(e, v, N).items():
                    close(record[name], value, f"{context}/{name}")
                both = np.isfinite(e) & np.isfinite(v)
                difference = v[both] - e[both]
                for name, event in zip(outcome_names[:4], [difference == 0, difference == 1,
                                                          difference == 2, difference >= 3]):
                    close(record[name], float(event.sum() / B), context + "/" + name)
                outcomes = record[outcome_names].to_numpy(float)
                require(np.all((outcomes >= 0) & (outcomes <= 1)) and abs(outcomes.sum() - 1) < 1e-10,
                        f"Outcome fractions do not partition the random orders: {context}")
                close(fit.ville_threshold, 1 / alpha, context + "/Ville_threshold")
                require(fit.log_eville_threshold <= -math.log(alpha) + 1e-8, f"eVille threshold exceeds Ville: {context}")
                signed = (record.theta_hat - fit.full_prior_location_raw_theta) / fit.scale
                close(record.conflict_signed_location_scale, signed, context + "/signed_conflict")
                close(record.conflict_abs_location_scale, abs(signed), context + "/absolute_conflict")
                v_mean = independent_summary(v, N)["mean_capped_at_N"]
                close(record.ville_mean_minus_flat, v_mean - flat_mean, context + "/ville_minus_flat")

    example_errors = []
    require(example_ts is not None, "Example target was not tested")
    for prior, frame in example.groupby("prior", sort=False):
        frame = frame.sort_values("n")
        fit = fits.loc[fits.target_lineage.eq(settings["example_target"]) & fits.prior.eq(prior)].iloc[0]
        pi = FlatPrior() if prior == FLAT else engine.Prior(fit.family, -float(fit.full_prior_location_raw_theta), float(fit.scale))
        ns = frame.n.to_numpy(int)
        require(np.array_equal(ns, np.arange(2, len(frame) + 2)), "Example prefix sizes are not consecutive")
        require(np.allclose(frame.t_raw_chronos, -example_ts, rtol=1e-10, atol=1e-10), "Example t statistics differ")
        refined = np.array([engine.log_mixture(float(t), int(n), pi, max(96, int(settings["quadrature_order"])))
                            for t, n in zip(example_ts, ns)])
        error = float(np.max(abs(refined - frame.log_evidence.to_numpy(float))))
        example_errors.append(error)
        require(error < 1e-7, f"Example evidence is unstable under quadrature refinement: {prior}")
        running = np.maximum.accumulate(frame.log_evidence.to_numpy(float))
        require(np.allclose(running, frame.log_running_evidence, rtol=1e-11, atol=1e-11), "Example running maxima differ")
        pe = frame.p_eville.to_numpy(float)
        require(np.isfinite(pe).all() and np.all((pe >= 0) & (pe <= 1)) and np.all(np.diff(pe) <= 1e-10),
                "Example eVille p-values invalid")
        if prior == FLAT:
            require(frame.p_ville.isna().all(), "Improper example has Ville p-values")
            require(np.allclose(pe, flat_calibration(running), rtol=1e-10, atol=1e-12), "Flat example calibration differs")
        else:
            require(np.allclose(frame.p_ville, np.exp(np.minimum(0, -running)), rtol=1e-10, atol=1e-12), "Example Ville p-values differ")

    core = validate_core()
    return {
        "passed": True,
        "tested_lineages": len(eligible),
        "main_lineages": len(main),
        "small_lineages": len(eligible - main),
        "paired_random_orders_per_lineage": B,
        "method_stopping_decisions_reproduced": checked_decisions,
        "method_summary_rows_reproduced": checked_summaries,
        "boundary_rows_checked": boundary_rows,
        "maximum_example_log_evidence_refinement_error": max(example_errors),
        "core_checks": core,
        "interpretation": "Numerical verification for these cases, not certified quadrature bounds or validation of Gaussianity, independence, or retrospective selection.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path(__file__).resolve().parent / "results")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--core-only", action="store_true")
    args = parser.parse_args()
    report = {"passed": True, "core_checks": validate_core()} if args.core_only else validate_results(args.results)
    path = args.output or args.results / "validation.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

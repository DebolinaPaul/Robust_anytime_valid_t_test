"""Shared sequential calculations and paired observation-order summaries.

Monitoring begins at n0=2. The computational response is X=-Chronos,
corresponding to the scientific alternative of negative mean Chronos.
Permutation 0 retains the supplied observation order; subsequent orders use
stable per-lineage seeds. The improper flat prior has eVille calibration only.
Nonrejection uses N+1 internally; restricted summaries cap stopping times at N.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import math

import numpy as np
from scipy.optimize import brentq

from engine import log_mixture, calibration_direct, direct_log_threshold
from improper_flat import FlatPrior


cached_log_threshold = lru_cache(maxsize=2048)(direct_log_threshold)


class DirectCalibration:
    """CDF-mixture calibration with cached values for repeated running maxima."""

    def __init__(self, prior, n0=2, quadrature_order=48):
        self.prior = prior
        self.n0 = n0
        self.quadrature_order = quadrature_order
        self._p_cache = {}
        self._threshold_cache = {}

    def pvalue(self, log_evidence):
        values = np.asarray(log_evidence, float)
        result = []
        for value in values.ravel():
            key = float(value)
            if key not in self._p_cache:
                self._p_cache[key] = calibration_direct(
                    key, self.prior, self.n0, self.quadrature_order)
            result.append(self._p_cache[key])
        result = np.asarray(result).reshape(values.shape)
        return float(result) if result.ndim == 0 else result

    def log_threshold(self, alpha=0.05):
        if alpha not in self._threshold_cache:
            self._threshold_cache[alpha] = direct_log_threshold(
                self.prior, alpha, self.n0, self.quadrature_order)
        return self._threshold_cache[alpha]

    def diagnostics(self, alpha=0.05):
        threshold = self.log_threshold(alpha)
        improper = isinstance(self.prior, FlatPrior)
        return {"n0": self.n0, "method": "analytic improper-flat calibration" if improper else "numerical noncentral-t CDF mixture",
                "quadrature_order": None if improper else self.quadrature_order,
                "log_eville_threshold": threshold, "eville_threshold": math.exp(threshold),
                "ville_threshold": None if improper else 1 / alpha,
                "positive_prior_mass": None if improper else self.prior.positive_mass}


def lineage_seed(seed, gene, sex, lineage):
    """Stable across platforms, fitting rules and subsets of requested folds."""
    key = f"{int(seed)}\0{gene}\0{sex}\0{lineage}".encode()
    return int.from_bytes(hashlib.sha256(key).digest()[:8], "little")


def permutation_t_statistics(x, permutations, seed):
    """Original and random orders, with stable prefix t statistics."""
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or x.size < 2 or not np.isfinite(x).all():
        raise ValueError("Provide at least two finite observations")

    if (
        not np.isfinite(permutations)
        or int(permutations) != permutations
        or permutations < 0
    ):
        raise ValueError("permutations must be a nonnegative integer")
    permutations = int(permutations)

    rng = np.random.default_rng(seed)
    N = x.size
    indices = np.vstack(
        [np.arange(N)]
        + [rng.permutation(N) for _ in range(permutations)]
    )
    observations = x[indices]
    ns = np.arange(2, N + 1)
    ts = np.empty((len(indices), N - 1), dtype=float)

    # Welford updates, vectorized across observation orders.
    origin = observations[:, 0]
    mean_offset = np.zeros(len(indices), dtype=float)
    m2 = np.zeros(len(indices), dtype=float)

    for j in range(1, N):
        n = j + 1
        offset = observations[:, j] - origin
        difference = offset - mean_offset
        mean_offset += difference / n
        m2 += difference * (offset - mean_offset)

        variance = m2 / (n - 1)
        if not np.isfinite(variance).all() or np.any(variance <= 0):
            raise ValueError(
                f"A prefix at n={n} has nonfinite or nonpositive variance"
            )

        ts[:, j - 1] = math.sqrt(n) * (
            (origin + mean_offset) / np.sqrt(variance)
        )

    if (
        not np.isfinite(ts).all()
        or np.max(np.abs(ts)) >= math.sinh(32)
    ):
        raise FloatingPointError(
            "Observed t statistics exceed the critical-boundary search range"
        )

    return indices, ts, ns


@lru_cache(maxsize=20000)
def critical_t(n, prior, log_threshold, quadrature_order=48):
    """Solve the monotone one-sided mixture boundary, including negative roots.

    asinh(t) gives a stable bracket from near-zero t to extreme tails. For a
    light-tailed prior, M_n(t) has a finite upper endpoint for fixed n; a
    threshold above that endpoint has no finite crossing. Heavy-tailed priors
    may have an infinite endpoint. An infinity return only means no crossing
    within the searched range, not a claim that a mathematical finite root
    never exists. Observed t values are required to lie strictly inside that
    range, so this sentinel gives exact crossing decisions for every evaluated
    ordering by monotonicity. Boundary CSVs label these sentinels explicitly.
    """
    def objective(z):
        return log_mixture(math.sinh(z), int(n), prior, quadrature_order) - log_threshold

    lower, upper = -16.0, 16.0
    flo, fhi = objective(lower), objective(upper)
    if not (np.isfinite(flo) and np.isfinite(fhi)):
        raise FloatingPointError(
            "Nonfinite mixture value while bracketing a critical t value.")
    while flo >= 0 and abs(lower) < 32:
        lower -= 4
        flo = objective(lower)
    while fhi < 0 and upper < 32:
        upper += 4
        fhi = objective(upper)
    if flo >= 0:
        return -math.inf
    if fhi < 0:
        return math.inf
    root = brentq(objective, lower, upper, xtol=1e-11, rtol=1e-12)
    return math.sinh(root)


def first_crossings(ts, boundaries, N):
    crossed = ts >= np.asarray(boundaries)
    # N+1 is only an internal censoring sentinel. Raw CSV tau is missing for NR.
    return np.where(crossed.any(axis=1), crossed.argmax(axis=1) + 2, N + 1)


def boundaries_for_evaluated_orders(ts, ns, prior, threshold_e, threshold_v, quadrature_order):
    """Compute only boundaries that can affect these paired stopping times.

    Once every fixed/random order has crossed a method's threshold, later
    boundaries cannot alter its first crossing. They are stored as NaN with
    an explicit not-needed status, rather than being called infinite or NR.
    This optimization does not change the full fixed-order evidence path.
    """
    e_boundaries = np.full(len(ns), np.nan)
    if threshold_v is None:
        # The improper flat comparator has no Ville computation or boundary array.
        e_crossed = np.zeros(len(ts), dtype=bool)
        for j, n in enumerate(ns):
            e_boundaries[j] = critical_t(
                n, prior, threshold_e, quadrature_order)
            e_crossed |= ts[:, j] >= e_boundaries[j]
            if e_crossed.all():
                break
        return e_boundaries, None
    v_boundaries = np.full(len(ns), np.nan)
    e_crossed = np.zeros(len(ts), dtype=bool)
    v_crossed = np.zeros(len(ts), dtype=bool)
    for j, n in enumerate(ns):
        v_boundaries[j] = critical_t(n, prior, threshold_v, quadrature_order)
        v_crossed |= ts[:, j] >= v_boundaries[j]
        if not e_crossed.all():
            e_boundaries[j] = (
                v_boundaries[j]
                if threshold_e == threshold_v
                else critical_t(n, prior, threshold_e, quadrature_order)
            )
            e_crossed |= ts[:, j] >= e_boundaries[j]
        if v_crossed.all():
            if not e_crossed.all():
                raise AssertionError(
                    "All orders crossed Ville but some failed to cross eVille.")
            break
    return e_boundaries, v_boundaries


def rejection_summary(times, N):
    times = np.asarray(times, dtype=int)
    events = times <= N
    observed = times[events]
    capped = np.minimum(times, N)
    k = int(math.ceil(0.5 * len(times)))
    empirical_median = int(np.sort(times)[k - 1])
    median_identified = empirical_median <= N
    return {
        "random_permutations": len(times),
        "rejection_count": int(events.sum()),
        "rejection_fraction": float(events.mean()),
        "mean_capped_at_N": float(capped.mean()),
        "median_capped_at_N": float(np.median(capped)),
        "mean_capped_at_N_over_N": float(capped.mean() / N),
        "median_capped_at_N_over_N": float(np.median(capped) / N),
        "unconditional_median_rejection_time": empirical_median if median_identified else None,
        "unconditional_median_status": "observed" if median_identified else "NR",
        "mean_first_rejection_time_if_all_orders_reject": float(times.mean()) if events.all() else None,
        "conditional_on_rejection_mean": float(observed.mean()) if len(observed) else None,
        "conditional_on_rejection_median": float(np.median(observed)) if len(observed) else None,
    }


def paired_summary(eville, ville, N):
    if ville is None:
        # Missing means the comparison was not performed, never a failed test.
        return {
            "comparison_status": "not_compared_eville_only",
            "random_permutations": len(eville),
            **{name: None for name in [
                "eville_earlier_fraction", "both_reject_fraction",
                "eville_only_rejects_fraction", "neither_rejects_fraction",
                "earlier_and_both_reject_fraction", "mean_restricted_observations_saved",
                "median_restricted_observations_saved", "mean_saved_conditional_on_both_rejecting",
                "median_saved_conditional_on_both_rejecting",
            ]},
        }
    e, v = np.asarray(eville), np.asarray(ville)
    if not np.all(e <= v):
        raise AssertionError(
            "eVille must reject no later than Ville for the same evidence path.")
    both = (e <= N) & (v <= N)
    difference = v[both] - e[both]
    capped_difference = np.minimum(v, N) - np.minimum(e, N)
    return {
        "comparison_status": "compared",
        "random_permutations": len(e),
        "eville_earlier_fraction": float(np.mean(e < v)),
        "both_reject_fraction": float(np.mean(both)),
        "eville_only_rejects_fraction": float(np.mean((e <= N) & (v > N))),
        "neither_rejects_fraction": float(np.mean((e > N) & (v > N))),
        "earlier_and_both_reject_fraction": float(np.mean((e < v) & both)),
        "mean_restricted_observations_saved": float(np.mean(capped_difference)),
        "median_restricted_observations_saved": float(np.median(capped_difference)),
        "mean_saved_conditional_on_both_rejecting": float(np.mean(difference)) if len(difference) else None,
        "median_saved_conditional_on_both_rejecting": float(np.median(difference)) if len(difference) else None,
    }

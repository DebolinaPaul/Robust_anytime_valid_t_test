"""Ordered one-sample t-tests and stored result summaries."""
from __future__ import annotations

from dataclasses import dataclass, field
from numbers import Real
import warnings

import numpy as np
import pandas as pd

from .priors import Prior, ev_prior
from .evidence import ev_log_mixture
from .calibration import Calibration, ev_calibration


def _finite_scalar(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real scalar")
    value = float(value)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be a finite real scalar")
    return value


def _observations(x, mu):
    try:
        values = np.asarray(x)
    except (TypeError, ValueError) as exc:
        raise ValueError("x must be a one-dimensional finite numeric sequence") from exc
    if values.ndim != 1 or values.size < 2 or values.dtype.kind not in "iuf":
        raise ValueError("x must be a one-dimensional finite numeric sequence of length at least two")
    values = values.astype(float, copy=True)
    if not np.all(np.isfinite(values)):
        raise ValueError("x must contain only finite observations; missing values are not removed")
    mu = _finite_scalar(mu, "mu")
    with np.errstate(over="ignore", invalid="ignore"):
        centered = values - mu
    if not np.all(np.isfinite(centered)):
        raise ValueError("Centering x at mu exceeds floating-point range")
    return centered


def _t_centered(centered):
    # Match the notebook: normalize before calculating the centered variance.
    # Recomputing each prefix avoids cancellation in a difference of squares.
    scale = np.max(np.abs(centered))
    if scale == 0:
        raise ValueError("Positive sample variance is required")
    z = centered / scale
    sd = np.std(z, ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("Positive sample variance is required")
    with np.errstate(over="ignore", invalid="ignore"):
        t = np.sqrt(len(z)) * np.mean(z) / sd
    if not np.isfinite(t):
        raise FloatingPointError("The t statistic exceeds floating-point range")
    return float(t)


def t_statistic(x, mu=0.0):
    """Student statistic for finite observations and null mean ``mu``.

    Computes ``sqrt(n) * mean(x - mu) / sd(x)`` with the sample standard
    deviation (``ddof=1``). Normalization protects against changes of units.
    """
    return _t_centered(_observations(x, mu))


def ev_log_evidence(x, prior=None, mu=0.0):
    """Log mixture evidence from observations; see ``ev_log_mixture`` for t input."""
    prior = ev_prior() if prior is None else prior
    if not isinstance(prior, Prior):
        raise TypeError("prior must be constructed with ev_prior()")
    y = _observations(x, mu)
    return float(ev_log_mixture(_t_centered(y), len(y), prior))


@dataclass
class TTestResult:
    """Stored one-sample t-test results; plotting does not rerun the test.

    ``None`` denotes no observed rejection in ``stopping_times``. For Lebesgue
    measure, Ville is unavailable; the summary distinguishes this case.
    ``path`` and ``to_dataframe()`` return copies of the stored trajectory.
    """

    prior: Prior
    calibration: Calibration
    mu: float
    alpha: float
    stopping_times: dict
    calibration_ambiguous_n: tuple
    relative_tolerance: float | None
    _path: pd.DataFrame = field(repr=False)

    @property
    def path(self):
        return self._path.copy(deep=True)

    def to_dataframe(self):
        """Return running p-values, evidence and numerical diagnostics."""
        return self.path

    def summary(self):
        """Return decisions and first rejection times as a pandas DataFrame.

        Monte Carlo SE is for the un-capped empirical calibration at the final
        running evidence, not a confidence interval or a validity certificate.
        Additional test settings and numerical diagnostics are in ``attrs``.
        """
        last = self._path.iloc[-1]
        mc = self.calibration.method == "monte_carlo"
        rows = []
        for method, col in (("eVille", "p_eville"), ("Ville", "p_ville")):
            available = method == "eVille" or self.prior.proper
            first = self.stopping_times[method]
            rows.append({
                "method": method,
                "final_running_p": float(last[col]) if available else np.nan,
                "first_rejection": first,
                "decision": ("Unavailable" if not available else
                             "No rejection observed" if first is None else "Rejected"),
                "calibration": (self.calibration.method if method == "eVille" else
                                "ville" if available else "unavailable"),
                "mc_samples": self.calibration.B if mc and method == "eVille" else 0,
                "final_mcse": float(last["mcse"]) if mc and method == "eVille" else np.nan,
            })
        result = pd.DataFrame(rows)
        result["first_rejection"] = pd.array(result["first_rejection"], dtype="Int64")
        result.attrs.update({
            "test": "one-sample t-test",
            "alternative": self.prior.alternative,
            "prior": repr(self.prior),
            "mu": self.mu,
            "alpha": self.alpha,
            "n": int(last["n"]),
            "n0": 2,
            "max_calibration_gap": float(self._path.calibration_gap.max()),
            "calibration_ambiguous_n": self.calibration_ambiguous_n,
            "monte_carlo_note": ("Finite Monte Carlo approximation; numerical brackets do not bound "
                                 "Monte Carlo sampling error." if mc else None),
        })
        return result

    def plot(self, method="both", ax=None, log_y=False, mark_rejections=True,
             title=None, legend=True):
        """Plot stored running p-values and return a matplotlib Axes.

        ``method`` is ``'both'``, ``'eVille'`` or ``'Ville'``. Ville cannot be
        plotted for Lebesgue measure. ``log_y`` uses a logarithmic vertical
        scale; exact floating-point zeros are displayed at the smallest
        positive normal float, without changing the stored results.
        """
        import matplotlib.pyplot as plt

        if method not in ("both", "eVille", "Ville"):
            raise ValueError("method must be 'both', 'eVille' or 'Ville'")
        if not self.prior.proper and method == "Ville":
            raise ValueError("Ville calibration is unavailable for Lebesgue measure")
        selected = ["eVille", "Ville"] if method == "both" else [method]
        if not self.prior.proper:
            selected = ["eVille"]
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 4.2), constrained_layout=True)
        mc = self.calibration.method == "monte_carlo"
        for name in selected:
            col = "p_eville" if name == "eVille" else "p_ville"
            values = self._path[col].to_numpy(copy=True)
            if log_y:
                values = np.maximum(values, np.finfo(float).tiny)
            colour = "#0072B2" if name == "eVille" else "#D55E00"
            label = "eVille (MC)" if name == "eVille" and mc else name
            ax.step(self._path.n, values, where="post", color=colour,
                    linestyle="-" if name == "eVille" else "--",
                    linewidth=1.8, label=label)
            if len(values) == 1:
                ax.scatter(self._path.n, values, color=colour, s=24)
            first = self.stopping_times[name]
            if mark_rejections and first is not None:
                idx = int(np.flatnonzero(self._path.n.to_numpy() == first)[0])
                ax.scatter([first], [values[idx]], color=colour, s=32, zorder=4)
        ax.axhline(self.alpha, color="0.4", linestyle=":", linewidth=1,
                   label=f"alpha = {self.alpha:g}")
        ax.set_xlabel("Number of observations")
        ax.set_ylabel("Running p-value (Monte Carlo eVille)" if mc else "Running anytime-valid p-value")
        ax.set_title(title if title is not None else
                     "One-sample t-test: Monte Carlo eVille calibration" if mc else
                     "Anytime-valid one-sample t-test")
        if log_y:
            ax.set_yscale("log")
        else:
            ax.set_ylim(-0.02, 1.02)
        if legend:
            ax.legend(frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        return ax

    def __repr__(self):
        note = (f"Monte Carlo calibration ({self.calibration.B:,} null draws; approximation)"
                if self.calibration.method == "monte_carlo" else
                f"{self.calibration.method} calibration")
        return (f"One-sample t-test | {self.prior.alternative} | alpha={self.alpha:g}\n"
                f"{self.prior!r}\n{note}\n" + self.summary().to_string(index=False))


def ev_ttest(x, prior=None, mu=0.0, alpha=0.05, n0=2, calibration=None,
             tolerance=1e-5, calibration_method="monte_carlo", null_t=None,
             relative_tolerance=1e-3):
    """Monitor an ordered one-sample t-test from the second observation.

    Assumes independent Gaussian observations with unknown common variance.
    The prior is on ``(mean - mu) / sigma`` and must be fixed independently of
    the monitored target data. One-sided support is chosen in ``ev_prior``.
    All prefixes must have positive sample variance. No observations are
    sorted, dropped or permuted. Only ``n0=2`` is currently implemented.

    The default uses the simulation notebook's empirical Monte Carlo
    calibration. Finite Monte Carlo calibration approximates the theoretical
    anytime-valid test; its numerical integration bounds do not control the
    sampling error of the null Monte Carlo sample. The deterministic option
    integrates the population calibration; Lebesgue calibration is analytic.

    A supplied calibration must match the exact prior and takes precedence
    over ``tolerance`` and ``calibration_method``. It can be reused across
    independent target datasets; ``null_t`` cannot also be supplied.
    """
    prior = ev_prior() if prior is None else prior
    if not isinstance(prior, Prior):
        raise TypeError("prior must be constructed with ev_prior()")
    alpha = _finite_scalar(alpha, "alpha")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be strictly between zero and one")
    if isinstance(n0, (bool, np.bool_)) or not isinstance(n0, Real) or n0 != 2:
        raise ValueError("Only n0=2 is implemented")
    if relative_tolerance is not None:
        relative_tolerance = _finite_scalar(relative_tolerance, "relative_tolerance")
        if not 0 < relative_tolerance < 1:
            raise ValueError("relative_tolerance must be None or strictly between zero and one")
    if calibration_method not in ("monte_carlo", "deterministic"):
        raise ValueError("calibration_method must be 'monte_carlo' or 'deterministic'")
    mu = _finite_scalar(mu, "mu")
    y = _observations(x, mu)
    n = np.arange(2, len(y) + 1)
    # Validate the complete observed path before doing expensive calibration.
    ts = []
    for k in n:
        try:
            ts.append(_t_centered(y[:k]))
        except (ValueError, FloatingPointError) as exc:
            raise type(exc)(f"At prefix n={k}: {exc}") from exc
    if calibration is None:
        calibration = ev_calibration(prior, tolerance=tolerance,
                                     method=calibration_method, null_t=null_t)
    elif null_t is not None:
        raise ValueError("Supply null_t to ev_calibration(), not alongside a calibration object")
    if not isinstance(calibration, Calibration) or calibration.prior != prior:
        raise ValueError("calibration must match the exact prior and alternative")
    log_evidence = np.asarray([ev_log_mixture(t, int(k), prior) for k, t in zip(n, ts)])
    if np.any(~np.isfinite(log_evidence)):
        raise FloatingPointError("Mixture log evidence is outside the supported numerical range")
    running = np.maximum.accumulate(log_evidence)
    cal = pd.DataFrame(calibration.evaluate_many(running, relative_tolerance=relative_tolerance))
    pe = np.minimum.accumulate(cal.p_eville.to_numpy())
    lower = np.minimum(pe, np.minimum.accumulate(cal.lower.to_numpy()))
    pv = np.exp(np.minimum(0.0, -running)) if prior.proper else np.full(len(n), np.nan)
    path = pd.DataFrame({"n": n, "log_evidence": log_evidence,
                         "running_log_evidence": running,
                         "p_eville": pe, "p_ville": pv,
                         "calibration_lower": lower,
                         "calibration_gap": np.maximum(0.0, pe - lower)})
    if calibration.method == "monte_carlo":
        for name in ("mc_estimate", "mc_lower", "mc_upper", "mc_gap", "mcse"):
            path[name] = cal[name].to_numpy()
        path["mc_samples"] = cal.B.to_numpy()
    ambiguous = tuple(int(k) for k in n[(lower <= alpha) & (pe > alpha)])
    if ambiguous:
        warnings.warn("The numerical calibration bracket straddles alpha at n="
                      + ", ".join(map(str, ambiguous[:8]))
                      + (", ..." if len(ambiguous) > 8 else "")
                      + ". Decisions use the numerical upper value. Tighten tolerance to "
                      "resolve numerical uncertainty; this does not remove Monte Carlo error.",
                      RuntimeWarning, stacklevel=2)
    stopping = {}
    for name, values in (("eVille", pe), ("Ville", pv)):
        hit = np.flatnonzero(values <= alpha)
        stopping[name] = int(n[hit[0]]) if len(hit) else None
    return TTestResult(prior=prior, calibration=calibration, mu=mu, alpha=alpha,
                       stopping_times=stopping, calibration_ambiguous_n=ambiguous,
                       relative_tolerance=relative_tolerance, _path=path)

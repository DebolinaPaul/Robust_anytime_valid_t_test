"""Null calibration for the Gaussian one-sample anytime-valid t-test.

The default Monte Carlo rule matches the accompanying simulation notebook.
Its brackets describe numerical error in a *finite empirical mean*, not
Monte Carlo error relative to the population calibration. The finite-sample
Monte Carlo approximation does not itself have an exact anytime guarantee.
"""

from __future__ import annotations

from functools import lru_cache
from importlib.resources import files
import gzip
import math
import numbers
import warnings

import numpy as np
import pandas as pd

from ._notebook_calibration import MCCalibration, NullSample, log_g_lebesgue
from .evidence import ev_log_mixture
from .priors import Prior, ev_prior


@lru_cache(maxsize=1)
def _reference_sample():
    resource = files("eville").joinpath("data/python_null_t_seed20260510.csv.gz")
    with resource.open("rb") as raw, gzip.open(raw, "rt") as stream:
        result = np.loadtxt(stream, skiprows=1, dtype=float)
    if result.shape != (100000,) or not np.isfinite(result).all():
        raise RuntimeError("The bundled null sample is incomplete or invalid")
    result.setflags(write=False)
    return result


def ev_null_sample(n=100000):
    """Return a copy of the first ``n`` shared null t statistics.

    The bundled 100,000 observations are generated from independent pairs of
    standard normals using NumPy's generator seed 20260510, exactly as in the
    paper notebook and R package. Returning this sample does not alter an RNG.
    Reuse it across priors; never construct calibration draws from target data.
    """
    if (isinstance(n, (bool, np.bool_)) or not isinstance(n, numbers.Integral)
            or not 2 <= n <= 100000):
        raise ValueError("n must be an integer between 2 and 100000")
    return _reference_sample()[:int(n)].copy()


def _real_vector(values, name, *, allow_infinite=False):
    raw = np.asarray(values)
    if (raw.ndim != 1 or raw.size == 0 or raw.dtype.kind not in "iuf"
            or np.isnan(raw).any()
            or (not allow_infinite and not np.isfinite(raw).all())):
        qualifier = "without NaN" if allow_infinite else "finite"
        raise ValueError(f"{name} must be a nonempty one-dimensional {qualifier} numeric array")
    return np.asarray(raw, dtype=float)


def _validate_relative(relative_tolerance):
    if relative_tolerance is not None and (
            isinstance(relative_tolerance, (bool, np.bool_))
            or not isinstance(relative_tolerance, numbers.Real)
            or not np.isfinite(relative_tolerance)
            or not 0 < relative_tolerance < 1):
        raise ValueError("relative_tolerance must be None or a number between zero and one")


def _unit_exp(value):
    return math.exp(min(0.0, value))


def _g_gaussian(log_c, scale):
    """Closed population calibration for a centered two-sided Gaussian."""
    b = 2.0 * scale * scale
    if not math.isfinite(b):
        raise ValueError("Gaussian scale is too large for analytic calibration")
    log_a = 0.5 * math.log1p(b)
    if log_c <= -log_a:
        return 1.0
    if log_c >= log_a:
        return math.exp(-log_c)
    q = min(1.0, max(0.0, math.expm1(log_a-log_c)/b))
    cosine, sine = math.sqrt(q), math.sqrt(max(0.0, 1.0-q))
    primitive = math.atan2(sine, math.exp(log_a)*cosine)
    return (2.0/math.pi)*(math.exp(-log_c)*primitive+math.asin(cosine))


def _affine_area(a, b, anchor_u, anchor_y, slope, floor=0.0):
    """Integrate the clipped affine function of sqrt(2)*sin(angle)."""
    cuts = [a, b]
    if slope != 0:
        for level in (floor, 1.0):
            root_u = anchor_u+(level-anchor_y)/slope
            if math.isfinite(root_u) and abs(root_u) < math.sqrt(2):
                root = math.asin(root_u/math.sqrt(2))
                if a < root < b:
                    cuts.append(root)
    cuts = sorted(set(cuts))
    total = 0.0
    for left, right in zip(cuts[:-1], cuts[1:]):
        width, mid = right-left, (right+left)/2
        u_mid = math.sqrt(2)*math.sin(mid)
        y_mid = anchor_y+slope*(u_mid-anchor_u)
        if y_mid <= floor:
            total += floor*width
            continue
        if y_mid >= 1:
            total += width
            continue
        q = width/2
        sinc_minus_one = (-q*q/6+q**4/120-q**6/5040 if abs(q) < 1e-3
                          else math.sin(q)/q-1)
        average = anchor_y+slope*((u_mid-anchor_u)+u_mid*sinc_minus_one)
        total += width*min(1.0, max(floor, average))
    return total


class _PopulationCalibration:
    """Convex upper/lower integration bounds for the population null law.

    Under the null T_2=tan(V), V uniform on (-pi/2, pi/2). M_2 is convex
    in U=sqrt(2)*sin(V). Chords bound it above on their own intervals;
    opposite-half extrapolated chords bound it below. Endpoint cells use
    the universal [0,1] bounds on the clipped integrand and are refined.
    """

    def __init__(self, prior, tolerance, max_intervals):
        self.prior, self.tolerance = prior, tolerance
        self.max_intervals = max_intervals
        self._values, self._cells = {}, None

    def _at(self, angle):
        if abs(angle) == math.pi/2:
            return math.inf
        if angle not in self._values:
            value = float(ev_log_mixture(math.tan(angle), 2, self.prior))
            if not math.isfinite(value):
                raise FloatingPointError("Nonfinite evidence in deterministic calibration")
            self._values[angle] = value
        return self._values[angle]

    def _cell(self, a, b):
        angles = (a, (a+b)/2, b)
        return angles, tuple(self._at(v) for v in angles)

    def _bounds(self, cell, log_c):
        angles, logs = cell
        widths = np.diff(angles)
        center_log = self._at(0.0)
        floors = []
        for i in range(2):
            # Convexity proves monotonicity beyond a nonnegative secant
            # from the central anchor (or a nonpositive secant to it).
            # This provides a rigorous constant lower bound in the tails.
            if angles[i] > 0 and math.isfinite(logs[i]) and logs[i] >= center_log:
                floor = _unit_exp(logs[i]-log_c)
            elif angles[i+1] < 0 and math.isfinite(logs[i+1]) and logs[i+1] >= center_log:
                floor = _unit_exp(logs[i+1]-log_c)
            else:
                floor = 0.0
            floors.append(floor)
        lower_rectangle = sum(w*f for w, f in zip(widths, floors))/math.pi
        if not all(math.isfinite(value) for value in logs):
            return np.array([lower_rectangle, (angles[-1]-angles[0])/math.pi])
        upper_rectangle = sum(width*_unit_exp(max(logs[i:i+2])-log_c)
                              for i, width in enumerate(widths))/math.pi
        u = math.sqrt(2)*np.sin(angles)
        if (max(logs)-log_c >= 600 or np.min(np.diff(u)) <= 1e-12
                or np.min(widths) <= 1e-7):
            return np.array([lower_rectangle, upper_rectangle])
        f = np.exp(np.asarray(logs)-log_c)
        slopes = np.diff(f)/np.diff(u)
        if not np.isfinite(slopes).all():
            return np.array([lower_rectangle, upper_rectangle])
        high = (_affine_area(angles[0], angles[1], u[0], f[0], slopes[0])
                + _affine_area(angles[1], angles[2], u[1], f[1], slopes[1]))/math.pi
        low = (_affine_area(angles[0], angles[1], u[1], f[1], slopes[1], floors[0])
               + _affine_area(angles[1], angles[2], u[1], f[1], slopes[0], floors[1]))/math.pi
        return np.array([max(lower_rectangle, low), min(upper_rectangle, high)])

    def evaluate(self, log_c):
        if not math.isfinite(log_c):
            p = float(log_c < 0)
            return p, p, 0
        if self._cells is None:
            edges = np.linspace(-math.pi/2, math.pi/2, 9)
            self._cells = [self._cell(a, b) for a, b in zip(edges[:-1], edges[1:])]
        cells = self._cells
        bounds = [self._bounds(cell, log_c) for cell in cells]
        ville = _unit_exp(-log_c)
        while True:
            total = np.sum(bounds, axis=0)
            lower, upper = float(total[0]), min(1.0, ville, float(total[1]))
            if lower > upper+5e-9:
                raise FloatingPointError("Population calibration bounds crossed; inspect evidence quadrature")
            lower = min(lower, upper)
            gap = upper-lower
            if gap <= self.tolerance or len(cells) >= self.max_intervals:
                break
            j = max(range(len(bounds)), key=lambda k: bounds[k][1]-bounds[k][0])
            (a, mid, b), _ = cells[j]
            if a == mid or mid == b:
                break
            left, right = self._cell(a, mid), self._cell(mid, b)
            cells[j], bounds[j] = left, self._bounds(left, log_c)
            cells.append(right)
            bounds.append(self._bounds(right, log_c))
        if gap > self.tolerance:
            warnings.warn(f"Population calibration gap {gap:.3g} exceeds tolerance "
                          f"{self.tolerance:.3g}; returning the conservative upper value. "
                          "Increase max_intervals to refine it.", RuntimeWarning, stacklevel=3)
        return upper, lower, len(cells)


class Calibration:
    """Reusable calibration for a fixed prior and monitoring from n=2.

    Construct with :func:`ev_calibration`. ``evaluate`` accepts a scalar
    natural-log evidence threshold, while ``evaluate_many`` accepts a vector.
    The Monte Carlo lower/upper columns bracket the empirical mean only;
    ``mcse`` is a pointwise descriptive Monte Carlo standard error.
    """

    def __init__(self, prior, tolerance, max_intervals, method, null_t):
        self.prior, self.n0 = prior, 2
        self.tolerance, self.max_intervals = tolerance, max_intervals
        self.B, self.null_t, self.null_sample_source = 0, None, None
        self._engine = None
        if not prior.proper:
            self.method, self.approximation = "analytic", False
        elif method == "monte_carlo":
            self.method, self.approximation = method, True
            sample = ev_null_sample() if null_t is None else _real_vector(null_t, "null_t")
            if len(sample) < 2:
                raise ValueError("null_t must contain at least two finite null t statistics")
            self.null_t = sample.copy()
            self.null_t.setflags(write=False)
            self.B = len(sample)
            self.null_sample_source = ("Python notebook: seed 20260510" if null_t is None
                                       else "supplied null_t")
            self._engine = MCCalibration(
                prior, NullSample(sample), tolerance=tolerance, max_intervals=max_intervals,
                log_evidence_fn=lambda t: ev_log_mixture(t, 2, prior))
        elif prior.family == "gaussian" and prior.location == 0 and prior.alternative == "two.sided":
            self.method, self.approximation = "analytic", False
        else:
            self.method, self.approximation = "deterministic", False
            self._engine = _PopulationCalibration(prior, tolerance, max_intervals)

    def evaluate(self, log_c, relative_tolerance=1e-3):
        """Calibrate one natural-log threshold and return numerical diagnostics."""
        _validate_relative(relative_tolerance)
        if (isinstance(log_c, (bool, np.bool_)) or not isinstance(log_c, numbers.Real)
                or math.isnan(log_c)):
            raise ValueError("log_c must be a real scalar without NaN")
        lc = float(log_c)
        ville = _unit_exp(-lc) if self.prior.proper else math.nan
        if self.method == "monte_carlo":
            r = self._engine.evaluate(lc, relative_tolerance=relative_tolerance)
            lower = min(r.lower, r.ville_bound)
            return dict(p_eville=r.used, lower=lower, gap=max(0.0, r.used-lower),
                        intervals=r.intervals, mc_estimate=r.estimate, mc_lower=r.lower,
                        mc_upper=r.upper, mc_gap=r.gap, mcse=r.mcse, B=r.B,
                        ville_bound=r.ville_bound)
        if not self.prior.proper:
            shift = 0.0 if self.prior.alternative == "two.sided" else math.log(4)
            p = math.exp(log_g_lebesgue(lc+shift))
            lower, intervals = p, 0
        elif self.method == "analytic":
            p = _g_gaussian(lc, self.prior.scale)
            lower, intervals = p, 0
        else:
            p, lower, intervals = self._engine.evaluate(lc)
        return dict(p_eville=p, lower=lower, gap=max(0.0, p-lower), intervals=intervals,
                    mc_estimate=math.nan, mc_lower=math.nan, mc_upper=math.nan,
                    mc_gap=math.nan, mcse=math.nan, B=0, ville_bound=ville)

    def evaluate_many(self, log_c, relative_tolerance=1e-3):
        """Return one diagnostic mapping per supplied log threshold, in input order."""
        values = _real_vector(log_c, "log_c", allow_infinite=True)
        _validate_relative(relative_tolerance)
        unique, inverse = np.unique(values, return_inverse=True)
        results = [self.evaluate(float(value), relative_tolerance) for value in unique]
        return [results[i].copy() for i in inverse]

    def log_threshold(self, alpha=0.05):
        """Return a log threshold whose implemented calibration is at most alpha.

        For Monte Carlo this inverts the numerical empirical upper calibration,
        not the exact population calibration. The returned threshold is finite.
        """
        if (isinstance(alpha, (bool, np.bool_)) or not isinstance(alpha, numbers.Real)
                or not math.isfinite(alpha) or not 0 < alpha < 1):
            raise ValueError("alpha must lie strictly between zero and one")
        if self.method == "monte_carlo":
            return self._engine.log_threshold(float(alpha))
        high, low, step = -math.log(alpha), -math.log(alpha)-2.0, 2.0
        for _ in range(80):
            if self.evaluate(high)["p_eville"] <= alpha:
                break
            high += step
            step *= 2
        else:
            raise FloatingPointError("Could not bracket the calibration threshold")
        step = 2.0
        for _ in range(80):
            if self.evaluate(low)["p_eville"] > alpha:
                break
            low -= step
            step *= 2
        else:
            raise FloatingPointError("Could not bracket the calibration threshold")
        for _ in range(100):
            mid = (low+high)/2
            if self.evaluate(mid)["p_eville"] <= alpha:
                high = mid
            else:
                low = mid
            if high-low <= 2e-7:
                return high
        raise FloatingPointError("Calibration threshold bisection did not converge")

    def __repr__(self):
        return f"Calibration(prior={self.prior!r}, method={self.method!r}, B={self.B})"


def ev_calibration(prior=None, tolerance=1e-5, max_intervals=4096,
                   method="monte_carlo", null_t=None):
    """Construct calibration matching a fixed prior and alternative.

    ``method='monte_carlo'`` uses the notebook's shared 100,000 null draws by
    default. ``method='deterministic'`` brackets the population null-angle
    integral, conditional on accurate evidence evaluations. Centered two-sided
    Gaussian deterministic calibration and all Lebesgue calibrations use
    analytic formulas. Neither numerical bracket is an interval-arithmetic
    certificate for every floating-point operation.
    """
    if prior is None:
        prior = ev_prior()
    if not isinstance(prior, Prior):
        raise TypeError("prior must be a Prior object returned by ev_prior")
    if (isinstance(tolerance, (bool, np.bool_)) or not isinstance(tolerance, numbers.Real)
            or not math.isfinite(tolerance) or not 0 < tolerance < 1):
        raise ValueError("tolerance must lie strictly between zero and one")
    if (isinstance(max_intervals, (bool, np.bool_))
            or not isinstance(max_intervals, numbers.Integral) or max_intervals < 8):
        raise ValueError("max_intervals must be an integer at least 8")
    if method not in ("monte_carlo", "deterministic"):
        raise ValueError("method must be 'monte_carlo' or 'deterministic'")
    if null_t is not None and (method != "monte_carlo" or not prior.proper):
        raise ValueError("null_t is only used for Monte Carlo calibration of a proper prior")
    return Calibration(prior, float(tolerance), int(max_intervals), method, null_t)


def ev_pvalues(log_evidence, calibration, relative_tolerance=None):
    """Evaluate calibration at log evidence and return a pandas DataFrame.

    Use cumulative maxima of log evidence for running p-values. Values at
    current evidence are pointwise calibrations. Monte Carlo diagnostics do
    not certify the approximation's population anytime-valid guarantee.
    """
    if not isinstance(calibration, Calibration):
        raise TypeError("calibration must be returned by ev_calibration")
    values = np.asarray(log_evidence)
    if values.ndim == 0:
        values = values.reshape(1)
    return pd.DataFrame(calibration.evaluate_many(values, relative_tolerance))

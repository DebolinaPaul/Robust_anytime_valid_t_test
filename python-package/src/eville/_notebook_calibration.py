"""Monte Carlo calibration using shared independent null t statistics.

Convexity provides numerical brackets for the actual empirical calibration
mean. These brackets control numerical integration error, not Monte Carlo
error. Finite Monte Carlo calibration is an approximation to the population
calibration and does not itself establish an exact anytime-valid guarantee.
"""

from dataclasses import dataclass
import math
import numpy as np
from scipy.optimize import brentq

from ._notebook_engine import log_mixture



def log_g_lebesgue(log_c):
    """Analytic log calibration for the two-sided full-line Lebesgue mixture."""
    lc = np.asarray(log_c,dtype=float)
    if np.isnan(lc).any():
        raise ValueError("log_c must not contain NaN")
    log_z = np.minimum(0.0,(0.5*math.log(math.pi)-lc)/2.0)
    result = np.empty_like(log_z)
    small = log_z < -20
    result[small] = math.log(4/math.pi)+log_z[small]
    z = np.exp(log_z[~small])
    result[~small] = np.log((2/math.pi)*(np.arcsin(z)+z*np.sqrt(np.maximum(0,1-z*z))))
    return float(result) if result.ndim == 0 else result


def log_threshold_lebesgue(alpha):
    """Invert the analytic Lebesgue calibration in log-threshold units."""
    if not 0 < alpha < 1:
        raise ValueError("Require 0<alpha<1")
    lower = 0.5*math.log(math.pi)
    upper = max(lower+2,2*math.log(4/math.pi/alpha)+lower+2)
    return float(brentq(lambda value: log_g_lebesgue(value)-math.log(alpha),
                        lower,upper,xtol=2e-12,rtol=2e-14))


class NullSample:
    """A reusable finite sample of null t statistics with one degree of freedom."""

    def __init__(self, t_null):
        t = np.asarray(t_null, dtype=float).reshape(-1)
        if t.size < 2 or not np.isfinite(t).all():
            raise ValueError("Provide at least two finite null t statistics")
        self.t = np.sort(t)
        self.B = int(t.size)
        extended_t = self.t.astype(np.longdouble)
        self.u = np.sqrt(np.longdouble(2))*(extended_t/np.hypot(1, extended_t))
        self.prefix_u = np.r_[np.longdouble(0), np.cumsum(self.u, dtype=np.longdouble)]
        self.prefix_u2 = np.r_[np.longdouble(0), np.cumsum(self.u*self.u, dtype=np.longdouble)]

    def centered_moments(self, start, stop, center):
        count = stop-start
        if count <= 0:
            return 0, np.longdouble(0), np.longdouble(0)
        if count <= 256 or self.u[stop-1]-self.u[start] < 1e-3:
            delta = self.u[start:stop]-center
            return count, delta.sum(dtype=np.longdouble), np.dot(delta, delta)
        first = self.prefix_u[stop]-self.prefix_u[start]
        second = self.prefix_u2[stop]-self.prefix_u2[start]
        return count, first-count*center, max(np.longdouble(0), second-2*center*first+count*center*center)

    def dkw_radius(self, delta=0.01):
        """Optional distribution-free MC error diagnostic, not a correction."""
        if not 0 < delta < 1:
            raise ValueError("delta must lie between zero and one")
        return math.sqrt(math.log(2/delta)/(2*self.B))


@dataclass(frozen=True)
class MCResult:
    estimate: float
    lower: float
    upper: float
    gap: float
    mcse: float
    used: float
    ville_bound: float
    intervals: int
    B: int


@dataclass(frozen=True)
class _SampleCell:
    indices: tuple
    logs: tuple


def _unit_exp(log_value):
    return math.exp(min(0.0, float(log_value)))


class MCCalibration:
    """Empirical extended-Ville calibration from a fixed independent sample.

    ``lower`` and ``upper`` bracket the finite-sample calibration mean.
    ``estimate`` and ``mcse`` use the pointwise midpoint of the numerical
    lower and upper functions. The MCSE is descriptive; it is not a
    simultaneous confidence interval or a confidence interval for a root.
    ``used`` is the numerical upper estimate capped by the exact ordinary
    Ville bound. This cap preserves comparison with Ville but does not turn
    the Monte Carlo calibration into an exact population calibration.
    """

    def __init__(self, prior, null_sample, tolerance=1e-6, max_intervals=4096,
                 log_evidence_fn=None, log_mean=0.0):
        if not 0 < tolerance < 1 or max_intervals < 8:
            raise ValueError("Require 0<tolerance<1 and max_intervals>=8")
        self.prior = prior
        self.sample = null_sample if isinstance(null_sample, NullSample) else NullSample(null_sample)
        self.tolerance = float(tolerance)
        self.max_intervals = int(max_intervals)
        self.log_mean = float(log_mean)
        self._callback = log_evidence_fn
        self._values = {}
        self._cells = None
        self._thresholds = {}
        self._tail_cache = {}
        self.max_gap = 0.0

    def _at(self, index):
        index = int(index)
        if index not in self._values:
            t = float(self.sample.t[index])
            value = (log_mixture(t, 2, self.prior) if self._callback is None
                     else self._callback(t))
            if not np.isfinite(value):
                raise FloatingPointError("Nonfinite evidence at a sampled null statistic")
            self._values[index] = float(value)
        return self._values[index]

    def _cell(self, start, stop):
        indices = (int(start), int((start+stop)//2), int(stop))
        return _SampleCell(indices, tuple(self._at(i) for i in indices))

    def _constant_stats(self, start, stop, lower, upper):
        count = stop-start
        midpoint = 0.5*(lower+upper)
        return np.array([count*lower, count*upper, count*midpoint*midpoint], dtype=np.longdouble)

    def _pair_stats(self, start, stop, anchor, low_value, low_slope,
                    high_value, high_slope, threshold):
        """Exact sample moments of two clipped affine functions on an interval."""
        if stop <= start:
            return np.zeros(3, dtype=np.longdouble)
        u = self.sample.u
        cuts = {start, stop}
        for value, slope in [(low_value, low_slope), (high_value, high_slope)]:
            if slope == 0:
                continue
            for level in (np.longdouble(0), threshold):
                difference = level-value
                if abs(difference) > 3*abs(slope):
                    continue
                root = anchor+difference/slope
                k = int(np.searchsorted(u, root, side='left'))
                if start < k < stop:
                    cuts.add(k)
        cuts = sorted(cuts)
        total = np.zeros(3, dtype=np.longdouble)
        for left, right in zip(cuts[:-1], cuts[1:]):
            center = u[(left+right-1)//2]

            def coefficients(value, slope):
                at_center = value+slope*(center-anchor)
                if at_center <= 0:
                    return np.longdouble(0), np.longdouble(0)
                if at_center >= threshold:
                    return np.longdouble(1), np.longdouble(0)
                return at_center/threshold, slope/threshold

            low_a, low_b = coefficients(low_value, low_slope)
            high_a, high_b = coefficients(high_value, high_slope)
            count, first, second = self.sample.centered_moments(left, right, center)
            lower = count*low_a+low_b*first
            upper = count*high_a+high_b*first
            mid_a, mid_b = (low_a+high_a)/2, (low_b+high_b)/2
            moment2 = count*mid_a*mid_a+2*mid_a*mid_b*first+mid_b*mid_b*second
            total += [min(count, max(0, lower)), min(count, max(0, upper)),
                      min(count, max(0, moment2))]
        return total

    def _cell_stats(self, cell, log_c):
        a, m, b = cell.indices
        logs = np.asarray(cell.logs, dtype=np.longdouble)
        u = self.sample.u[[a,m,b]]
        if u[0] == u[2] and self.sample.t[a] != self.sample.t[b]:
            # Near the invariant boundary, distinct finite t values may
            # round to the same U on platforms where longdouble is float64.
            # Evaluate the actual observations rather than declaring a tie.
            values = np.asarray([_unit_exp(self._at(i)-log_c)
                                 for i in range(a,b)],dtype=np.longdouble)
            return np.array([values.sum(),values.sum(),np.dot(values,values)],
                            dtype=np.longdouble)
        if b-a == 1 or u[0] == u[2]:
            value = _unit_exp(logs[0]-log_c)
            return self._constant_stats(a,b,value,value)
        if u[0] == u[1] or u[1] == u[2]:
            return self._constant_stats(a,b,0.0,_unit_exp(max(logs)-log_c))
        shift = max(float(np.max(logs)), float(log_c))
        f = np.exp(logs-np.longdouble(shift))
        threshold = np.exp(np.longdouble(log_c)-shift)
        slopes = np.diff(f)/np.diff(u)
        if not np.isfinite(slopes).all():
            raise FloatingPointError("Nonfinite affine slope in Monte Carlo calibration")
        if threshold < np.longdouble('1e-250'):
            # Resolve extreme relative clipping by whole-interval bounds.
            left_floor = min(f[1], f[1]+slopes[1]*(u[0]-u[1]))
            right_floor = min(f[1], f[1]+slopes[0]*(u[2]-u[1]))
            left_lower = float(left_floor > 0 and np.log(left_floor)+shift >= log_c)
            right_lower = float(right_floor > 0 and np.log(right_floor)+shift >= log_c)
            left_upper = _unit_exp(max(logs[0],logs[1])-log_c)
            right_upper = _unit_exp(max(logs[1],logs[2])-log_c)
            return (self._constant_stats(a,m,left_lower,left_upper)
                    +self._constant_stats(m,b,right_lower,right_upper))
        # Secants bound a convex function above on their own halves;
        # extrapolating the opposite secant gives the lower support.
        left = self._pair_stats(a,m,u[1],f[1],slopes[1],f[1],slopes[0],threshold)
        right = self._pair_stats(m,b,u[1],f[1],slopes[0],f[1],slopes[1],threshold)
        return left+right

    def evaluate(self, log_c, relative_tolerance=None, _absolute_tolerance=None):
        if relative_tolerance is not None and not 0 < relative_tolerance < 1:
            raise ValueError("relative_tolerance must lie between zero and one")
        lc = float(log_c)
        absolute_tolerance = (self.tolerance if _absolute_tolerance is None
                              else float(_absolute_tolerance))
        if math.isnan(lc):
            raise ValueError("log_c must not be NaN")
        if not math.isfinite(lc):
            value = 1.0 if lc < 0 else 0.0
            return MCResult(value,value,value,0.0,0.0,value,value,0,self.sample.B)
        endpoint_max = max(self._at(0),self._at(self.sample.B-1))
        if lc > endpoint_max:
            factor = math.exp(endpoint_max-lc)
            anchor_tolerance = math.exp(min(0.0,math.log(absolute_tolerance)+lc-endpoint_max))
            cached = self._tail_cache.get(relative_tolerance)
            if cached is None or cached.gap > anchor_tolerance:
                self._tail_cache[relative_tolerance] = self.evaluate(
                    endpoint_max,relative_tolerance=relative_tolerance,
                    _absolute_tolerance=anchor_tolerance)
            anchor = self._tail_cache[relative_tolerance]
            lower,upper = anchor.lower*factor,anchor.upper*factor
            ville = _unit_exp(self.log_mean-lc)
            if _absolute_tolerance is None:
                self.max_gap = max(self.max_gap,anchor.gap*factor)
            return MCResult(anchor.estimate*factor,lower,upper,anchor.gap*factor,
                            anchor.mcse*factor,min(upper,ville),ville,
                            anchor.intervals,self.sample.B)
        if self._cells is None:
            edges = np.unique(np.linspace(0,self.sample.B-1,9,dtype=int))
            self._cells = [self._cell(a,b) for a,b in zip(edges[:-1],edges[1:])]
        cells = self._cells
        stats = [self._cell_stats(cell,lc) for cell in cells]
        # Each cell includes its left endpoint and excludes its right one.
        # Include the single maximum observation separately.
        last = np.longdouble(_unit_exp(self._at(self.sample.B-1)-lc))
        total = np.sum(stats,axis=0,dtype=np.longdouble)+[last,last,last*last]
        B = self.sample.B
        while True:
            lower,upper = float(total[0]/B),float(total[1]/B)
            if lower > upper+1e-9:
                raise FloatingPointError("Empirical calibration bounds crossed")
            lower = min(lower,upper)
            gap = max(0.0,upper-lower)
            requested = absolute_tolerance
            if relative_tolerance is not None and upper > 0:
                requested = min(requested,relative_tolerance*upper)
            if gap <= requested:
                break
            if len(cells) >= self.max_intervals:
                raise FloatingPointError(f"Empirical calibration gap {gap:.3g} exceeds {requested:.3g}; increase max_intervals")
            j = max(range(len(stats)),key=lambda k: stats[k][1]-stats[k][0])
            a,m,b = cells[j].indices
            if m == a or m == b:
                raise FloatingPointError("Sample interval cannot be refined further")
            left,right = self._cell(a,m),self._cell(m,b)
            ls,rs = self._cell_stats(left,lc),self._cell_stats(right,lc)
            total += ls+rs-stats[j]
            cells[j],stats[j] = left,ls
            cells.append(right)
            stats.append(rs)
        estimate = 0.5*(lower+upper)
        second = float(total[2]/B)
        mcse = math.sqrt(max(0.0,second-estimate*estimate)/(B-1))
        ville = _unit_exp(self.log_mean-lc)
        if _absolute_tolerance is None:
            self.max_gap = max(self.max_gap,gap)
        return MCResult(estimate,lower,upper,gap,mcse,min(upper,ville),ville,len(cells),B)

    def log_upper(self, log_c, relative_tolerance=None):
        """Log of the numerical upper MC calibration capped by Ville."""
        lc = float(log_c)
        endpoint_max = max(self._at(0),self._at(self.sample.B-1))
        if math.isfinite(lc) and lc > endpoint_max:
            # All sampled evidence is below c. The empirical calibration
            # equals its value at the endpoint maximum times exp(anchor-lc).
            self.evaluate(lc,relative_tolerance=relative_tolerance)
            anchor = self._tail_cache[relative_tolerance]
            return min(math.log(anchor.upper)+endpoint_max-lc,self.log_mean-lc,0.0)
        value = self.evaluate(lc,relative_tolerance=relative_tolerance).used
        return math.log(value) if value > 0 else -math.inf

    def log_threshold(self, alpha, log_tolerance=2e-7):
        if not 0 < alpha < 1:
            raise ValueError("Require 0<alpha<1")
        alpha = float(alpha)
        if alpha in self._thresholds:
            return self._thresholds[alpha]
        high = float(np.nextafter(self.log_mean-math.log(alpha),math.inf))
        low,step = high-2.0,2.0
        for _ in range(80):
            if self.evaluate(low).used > alpha:
                break
            high = low
            step *= 2
            low -= step
        else:
            raise FloatingPointError("Could not bracket empirical log threshold")
        for _ in range(80):
            middle = (low+high)/2
            if self.evaluate(middle).used <= alpha:
                high = middle
            else:
                low = middle
            if high-low <= log_tolerance:
                break
        if self.evaluate(high).used > alpha*(1+1e-12):
            raise FloatingPointError("Returned empirical threshold failed its check")
        self._thresholds[alpha] = high
        return high

    threshold = log_threshold

    def evaluate_many(self, log_c, relative_tolerance=None):
        values = np.asarray(log_c,dtype=float)
        unique,inverse = np.unique(values,return_inverse=True)
        results = [self.evaluate(x,relative_tolerance=relative_tolerance) for x in unique]
        return [results[i] for i in inverse]


Calibrator = MCCalibration

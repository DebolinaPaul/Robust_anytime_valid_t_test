"""Fixed-sample normal-model confidence intervals for mu / sigma.

Under iid N(mu, sigma**2) sampling, T=sqrt(n)*mean/s follows a
noncentral t distribution with n-1 degrees of freedom and noncentrality
sqrt(n)*mu/sigma. Inverting that distribution gives an equal-tailed
confidence interval. These are NOT confidence sequences and do not have
simultaneous coverage across lineages. "Exact" refers to the sampling
distribution, not to exact arithmetic: SciPy evaluates the CDF and a
bracketed numerical solver inverts it. Failure raises an exception; there
is no normal-approximation fallback.
"""
from __future__ import annotations

import math
import warnings

from scipy.optimize import brentq
from scipy.stats import nct, norm


def standardized_mean_ci(t_stat, n, confidence=0.95):
    """Return the fixed-sample equal-tailed CI for theta=mu/sigma.

    Parameters
    ----------
    t_stat : float
        sqrt(n) * sample_mean / sample_standard_deviation (ddof=1).
    n : int
        Sample size, at least two.
    confidence : float
        Confidence level strictly between zero and one.

    Negative statistics are handled by sign reflection. This avoids the
    inaccurate negative-t/negative-noncentrality regime encountered in
    some noncentral-t CDF implementations.
    """
    if isinstance(n, bool) or not math.isfinite(float(n)) or int(n) != n or n < 2:
        raise ValueError("n must be a finite integer at least 2")
    t_stat, confidence = float(t_stat), float(confidence)
    if not math.isfinite(t_stat):
        raise ValueError("t_stat must be finite")
    if not math.isfinite(confidence) or not 0.0 < confidence < 1.0:
        raise ValueError("confidence must lie strictly between 0 and 1")
    if t_stat < 0:
        low, high = standardized_mean_ci(-t_stat, n, confidence)
        return -high, -low

    df = int(n) - 1
    alpha = 1.0 - confidence

    def cdf(noncentrality):
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            value = float(nct.cdf(t_stat, df, noncentrality))
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise ArithmeticError("Noncentral-t CDF failed during CI inversion")
        return value

    def invert(probability):
        # F(t; lambda) is continuous and strictly decreasing in lambda.
        # First locate the root's sign, avoiding unnecessarily large
        # negative noncentralities (where some SciPy CDFs return NaN).
        # For t>=0, F(t; lambda)>=Phi(-lambda); consequently
        # -Phi^{-1}(p)-1 is a strict lower bracket whenever the root is
        # negative. The positive initial bracket is expanded as needed.
        at_zero = cdf(0.0)
        if at_zero < probability:
            left, right = -float(norm.ppf(probability)) - 1.0, 0.0
        else:
            left = 0.0
            right = max(8.0, t_stat * (1.0 + math.sqrt(2.0 / df)) + 8.0)
        for _ in range(100):
            if not math.isfinite(left) or not math.isfinite(right):
                raise ArithmeticError("Noncentrality bracket overflow")
            if cdf(left) >= probability >= cdf(right):
                break
            if at_zero < probability:
                raise ArithmeticError("The exact negative-root bracket failed numerically")
            right *= 2.0
        else:
            raise ArithmeticError("Could not bracket noncentral-t CI endpoint")
        root = float(brentq(lambda lam: cdf(lam) - probability,
                            left, right, xtol=1e-11, rtol=1e-12,
                            maxiter=200))
        residual = abs(cdf(root) - probability)
        if residual > 1e-8:
            raise ArithmeticError(f"CI endpoint CDF residual is {residual:g}")
        return root / math.sqrt(n)

    lower = invert(1.0 - alpha / 2.0)
    upper = invert(alpha / 2.0)
    if lower > upper:
        raise ArithmeticError("Noncentral-t inversion returned reversed endpoints")
    return lower, upper

"""Exact half-flat standardized-effect mixture and n0=2 extended Ville.

For X_i ~ N(mu,sigma^2), delta=mu/sigma, H0 delta<=0 vs delta>0,
use the sigma-finite measure d(delta) on (0,infinity). This measure is NOT
flat in unstandardized mu, is NOT a probability distribution, and has no
ordinary Ville calibration with unit mean. Starting monitoring at n0=2,
its clipped-mean extended-Ville calibration is exact up to floating point.

Writing nu=n-1 and V~chi-square(nu), integration over delta gives
  integral f_nct(t;nu,sqrt(n)*delta) ddelta
    = E[sqrt(V/nu)*Phi(t*sqrt(V/nu))]/sqrt(n).
Tilting the chi-square density by sqrt(V) and dividing by the central t
pdf yields
  M_n(t) = sqrt(2*pi/n)*(1+t*t/(n-1))**(n/2)
           * F_t[n](t*sqrt(n/(n-1))).
At n=2 this reduces to a*(1+exp(2*asinh(t))), a=sqrt(pi)/4.
Since T_2 is standard Cauchy under the boundary null,
  g(c)=E[min(1,M_2/c)] = (2*phi+sin(2*phi))/pi,
  phi=arcsin(sqrt(a/c)) for c>a; g(c)=1 for c<=a.

Conditional Tonelli gives the extended martingale property; E_0[M_n]=
int_0^infinity E_0[L_n(delta)] ddelta=infinity. Extended Ville remains
valid: set B={M_2<c}, so 1_B*M_n is an integrable martingale for n>=2.
Ordinary Ville on this truncated-initial-event process, plus P(B^c), gives
P(sup_{n>=2} M_n>=c)<=E[min(1,M_2/c)]. This is Theorem 4.1 of Wang and
Ramdas, arXiv:2304.01163; their Lemma 5.4 establishes sigma-finite mixing.
The scale-invariant Gaussian likelihood ratios are Corollary 4.6 of
Wang and Ramdas, arXiv:2310.03722. The one-sided closed form above is
obtained by the displayed integration, not assumed from a diffuse limit.
"""
from __future__ import annotations
import math
from dataclasses import dataclass, field
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from scipy.special import betaln, betainc, hyp2f1
from scipy.stats import t as student_t

FLAT_LABEL = 'Improper flat'
_LOG_A = .5*math.log(math.pi)-math.log(4.)


@dataclass(frozen=True)
class FlatPrior:
    """Lebesgue measure on positive standardized effects; never normalize it.

    scale=1 records the arbitrary unit thickness, not a probability scale.
    Multiplying thickness and evidence by k leaves calibrated p-values and
    rejection events unchanged when the calibration is scaled accordingly.
    """
    family: str = field(default=FLAT_LABEL, init=False)
    location: float = field(default=0., init=False)
    scale: float = field(default=1., init=False)
    proper: bool = field(default=False, init=False)
    positive_mass: float = field(default=math.inf, init=False)
    log_positive_mass: float = field(default=math.inf, init=False)

    def logpdf(self, delta):
        return np.where(np.asarray(delta, dtype=float) > 0, 0., -np.inf)

    def untruncated_logpdf(self, delta):
        return self.logpdf(delta)

    def pdf(self, delta):
        return np.exp(self.logpdf(delta))


def log_flat_mixture(t, n):
    """log ∫_0^∞ [f_nct(t;n−1,sqrt(n)*delta)/f_t(t;n−1)] d(delta).

    Positive-tail formula uses a Student-t CDF. Negative-tail formula uses
    its equivalent hypergeometric expression to cancel divergent powers
    analytically before numerical evaluation.
    """
    if int(n) != n or n < 2:
        raise ValueError('n must be an integer >= 2')
    n = int(n)
    ts = np.asarray(t, dtype=float)
    if not np.isfinite(ts).all():
        raise ValueError('t must be finite')
    if n == 2:
        result = _LOG_A+np.logaddexp(0., 2*np.arcsinh(ts))
    else:
        scaled = ts/math.sqrt(n-1)
        with np.errstate(over='ignore'):
            logz = -np.log1p(scaled*scaled)
        huge = ~np.isfinite(logz)
        logz = np.where(huge, math.log(n-1)-2 *
                        np.log(np.hypot(math.sqrt(n-1), ts)), logz)
        z = np.minimum(1., np.exp(logz))
        result = np.empty_like(ts)
        pos = ts >= 0
        result[pos] = (.5*math.log(2*math.pi/n)-.5*n*logz[pos]
                       + np.log1p(-.5*betainc(n/2, .5, z[pos])))
        zz = z[~pos]
        hh = hyp2f1(n/2, .5, n/2+1, zz)
        for i in np.flatnonzero(~np.isfinite(hh)):
            zval = float(zz.flat[i])
            a = n/2
            hh.flat[i] = quad(lambda v: 2*v*math.exp(-v*v)/math.sqrt((1-zval)+zval *
                              (-math.expm1(-v*v/a))), 0., np.inf, epsabs=1e-10, epsrel=1e-12)[0]
        result[~pos] = (.5*math.log(2*math.pi/n)-math.log(n)-betaln(n/2, .5)
                        + np.log(hh))
        # Near zero, evaluating 1-z through incomplete-beta arguments loses
        # precision. Direct t.logcdf is stable here and avoids that loss.
        near = np.abs(ts) <= 10.
        result[near] = (.5*math.log(2*math.pi/n)-.5*n*logz[near]
                        + student_t.logcdf(ts[near]*math.sqrt(n/(n-1)), df=n))
    if not np.isfinite(result).all():
        raise FloatingPointError('The analytic flat mixture was nonfinite')
    return float(result) if result.ndim == 0 else result


def flat_calibration(log_c, n0=2):
    """g(c)=E_0[min(1,M_2/c)], accepting log(c)."""
    if n0 != 2:
        raise ValueError('This analytic calibration is derived only for n0=2')
    values = np.asarray(log_c, dtype=float)
    if np.isnan(values).any():
        raise ValueError('log_c must not be NaN')
    s = np.exp(np.minimum(0., .5*(_LOG_A-values)))
    phi = np.arcsin(s)
    out = (2*phi+np.sin(2*phi))/math.pi
    return float(out) if out.ndim == 0 else out


def flat_log_threshold(alpha=.05, n0=2):
    """Log threshold using a root whose scale does not vanish with alpha."""
    alpha = float(alpha)
    if n0 != 2 or not math.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("Require n0=2 and 0 < alpha < 1")

    # Set r = phi/alpha. The root lies in [pi/4, pi/2].
    # NumPy defines sinc(x) = sin(pi*x)/(pi*x).
    def objective(r):
        return (2 * r / math.pi) * (
            1 + np.sinc(2 * alpha * r / math.pi)
        ) - 1

    r = brentq(
        objective, math.pi / 4, math.pi / 2,
        xtol=1e-14, rtol=1e-14
    )

    log_sin_phi = (
        math.log(alpha)
        + math.log(r)
        + math.log(float(np.sinc(alpha * r / math.pi)))
    )
    return _LOG_A - 2 * log_sin_phi


def flat_running_pvalues(t_values, n_values):
    """Flat-prior evidence and running p-values, monitoring from n=2."""
    t_values = np.asarray(t_values, dtype=float)
    supplied_n = np.asarray(n_values)  # Validate before integer conversion.

    if (t_values.ndim != 1
            or supplied_n.shape != t_values.shape
            or t_values.size == 0
            or not np.isfinite(t_values).all()):
        raise ValueError(
            "Provide equal nonempty vectors and finite t values"
        )

    ns = np.arange(2, 2 + t_values.size, dtype=int)
    if not np.array_equal(supplied_n, ns):
        raise ValueError("n_values must be exactly 2, 3, ..., N")

    log_m = np.array([
        log_flat_mixture(t, int(n))
        for t, n in zip(t_values, ns)
    ])
    return log_m, flat_calibration(np.maximum.accumulate(log_m))

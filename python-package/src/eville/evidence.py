"""Invariant likelihood ratios and proper/improper effect-size mixtures.

The two-sided proper-prior engine is the simulation notebook's engine.
Directional mixtures use effect-space quadrature with refinement and
an explicit likelihood-envelope bound on the omitted tail, as in R.
"""
from functools import lru_cache
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import betaln, logsumexp
from scipy import stats
from .priors import Prior, ev_prior, _logpdf
from . import _notebook_engine as _engine

_NAMES = {'gaussian': 'Gaussian', 'cauchy': 'Cauchy', 'student_t5': 'Student-t5',
          'laplace': 'Laplace', 'horseshoe': 'Horseshoe'}


def _half_grid(upper, u, n, d, location, scale):
    like_scale = 1/math.sqrt(d if u > 0 else n)
    like_center = u*math.sqrt((n-1)/(n*d)) if u > 0 else 0.0
    centers = np.unique([0.0, max(location, 0.0), like_center])
    scales = np.unique([scale, like_scale])
    offsets = np.r_[-2.0**np.arange(6, -5, -1), 0.0, 2.0**np.arange(-4, 7)]
    points = np.concatenate([center+s*offsets for center in centers for s in scales])
    low = min(scales)/16
    if 0 < low < upper:
        log_ratio = math.log(upper)-math.log(low)
        count = max(2, math.ceil(2*log_ratio/math.log(2))+1)
        if count > 4096:
            raise FloatingPointError('The prior/likelihood scale separation is too large.')
        points = np.r_[points, np.exp(np.linspace(math.log(low), math.log(upper), count))]
    points = points[np.isfinite(points) & (points > 0) & (points < upper)]
    return np.unique(np.r_[0.0, upper, points])


def _half_once(t, upper, u, n, d, family, location, scale, log_mass, order):
    breaks = _half_grid(upper, u, n, d, location, scale)
    left, right = breaks[:-1], breaks[1:]
    half = (right-left)/2
    mid = left+half
    q, log_w = _engine._legendre(order)
    theta = mid[:, None]+half[:, None]*q
    lw = np.log(half[:, None])+np.broadcast_to(log_w, theta.shape)
    lp = _logpdf(theta, family, location, scale)
    if family == 'horseshoe':
        # Use theta = location +/- r**4 on *every* interval on each side.
        # Only transforming the two adjacent intervals can leave a nearby
        # grid interval with a nearly singular endpoint and slow convergence.
        # When location is inside the integration range, the grid contains it,
        # so no interval crosses the spike.
        for j in range(len(left)):
            direction = 1 if left[j] >= location else -1
            distances = sorted((abs(left[j]-location), abs(right[j]-location)))
            r_left, r_right = np.sqrt(np.sqrt(distances))
            r_half = (r_right-r_left)/2
            if r_half <= 0:
                # An ordinary tiny interval far from the spike is smooth;
                # retain its original mapping if transformed endpoints round.
                continue
            r = r_left+r_half*(1+q)
            offset = direction*r**4
            theta[j] = location+offset
            lw[j] = log_w+math.log(r_half)+math.log(4)+3*np.log(r)
            lp[j] = _logpdf(offset, family, 0.0, scale)
    lr = _engine.log_likelihood_ratio(t, n, theta)
    return float(logsumexp(lr+lp-log_mass+lw))


def _half_mixture(t, n, prior, direction, rtol=1e-8):
    t = direction*t
    u, log_d = _engine._invariant(t, n)
    d = math.exp(log_d)
    if d <= 0:
        raise FloatingPointError('Residual energy is outside numerical range.')
    location = direction*prior.location
    log_mass = prior.log_positive_mass if direction == 1 else prior.log_negative_mass
    log_c = 0.0 if u <= 0 else n/2*math.log(2*n/d)
    kappa = n/2 if u <= 0 else n*d/(2*(2*n-d))
    upper = max(1.0, 4/math.sqrt(d if u > 0 else n))
    for _ in range(60):
        coarse = _half_once(t, upper, u, n, d, prior.family, location,
                            prior.scale, log_mass, 32)
        fine = _half_once(t, upper, u, n, d, prior.family, location,
                          prior.scale, log_mass, 64)
        discrepancy = abs(np.expm1(coarse-fine))
        if not math.isfinite(fine):
            raise FloatingPointError('The directional mixture is outside numerical range.')
        for order in (128, 256, 512):
            if math.isfinite(discrepancy) and discrepancy <= rtol/4:
                break
            refined = _half_once(t, upper, u, n, d, prior.family, location,
                                 prior.scale, log_mass, order)
            discrepancy = abs(np.expm1(fine-refined))
            fine = refined
        if not math.isfinite(discrepancy) or discrepancy > rtol/4:
            raise FloatingPointError('Directional mixture quadrature did not converge under order refinement.')
        log_tail = log_c-kappa*upper**2
        if log_tail-fine <= math.log(rtol/4):
            return fine
        required = math.sqrt(max(0.0, (log_c-fine-math.log(rtol/8))/kappa))
        upper = max(1.5*upper, 1.05*required)
        if not math.isfinite(upper):
            break
    raise FloatingPointError('Directional mixture omitted-tail bound did not converge.')


def _log_half_lebesgue(t, n, direction):
    t = direction*t
    if n == 2:
        return 0.5*math.log(math.pi)-math.log(4)+float(np.logaddexp(0, 2*math.asinh(t)))
    if t >= -10:
        return _engine.log_lebesgue(t, n)+float(stats.t.logcdf(t*math.sqrt(n/(n-1)), df=n))
    _, log_d = _engine._invariant(t, n)
    z = math.exp(log_d)/n
    a = n/2
    # Cancel opposing large powers before evaluating the negative tail.
    def integrand(q):
        return math.exp(-q)/math.sqrt((1-z)+z*(-math.expm1(-q/a)))
    value, error = quad(integrand, 0, np.inf, epsrel=1e-10, epsabs=0, limit=500)
    if value <= 0 or not math.isfinite(value) or error > 1e-9*value:
        raise FloatingPointError('Half-Lebesgue quadrature failed its error check.')
    return 0.5*math.log(2*math.pi/n)-math.log(n)-betaln(n/2, 0.5)+math.log(value)


@lru_cache(maxsize=50000)
def _scalar_log_mixture(t, n, prior):
    _engine._invariant(t, n)
    if not prior.proper:
        if prior.alternative == 'two.sided':
            return _engine.log_lebesgue(t, n)
        return _log_half_lebesgue(t, n, 1 if prior.alternative == 'greater' else -1)
    if prior.alternative == 'two.sided':
        p = _engine.Prior(_NAMES[prior.family], prior.location, prior.scale)
        return _engine.log_mixture(t, n, p)
    return _half_mixture(t, n, prior, 1 if prior.alternative == 'greater' else -1)


def ev_log_mixture(t, n, prior=None):
    """Log mixture evidence at a t statistic and sample size.

    Parameters
    ----------
    t : float or array-like
        Finite one-sample t statistic(s), using the sample SD with ddof=1.
    n : int
        Common sample size, at least two; positive sample variance is required.
    prior : Prior, optional
        Fixed standardized-effect prior, including one-sided or Lebesgue.

    For raw observations use ``ev_ttest``. A statistic's density ratio is
    a martingale here because it is the invariant likelihood ratio for
    the Gaussian model, not because arbitrary t-density ratios are martingales.
    """
    prior = ev_prior() if prior is None else prior
    if not isinstance(prior, Prior):
        raise TypeError('prior must be constructed with ev_prior().')
    n = _engine._check_n(n)
    values = np.asarray(t, dtype=float)
    if np.any(~np.isfinite(values)):
        raise ValueError('t must be finite, with positive sample variance.')
    answer = np.array([_scalar_log_mixture(float(v), n, prior) for v in values.ravel()]).reshape(values.shape)
    return float(answer) if answer.ndim == 0 else answer


def log_likelihood_ratio(t, n, theta):
    """Log invariant likelihood ratio at signed standardized effects theta."""
    return _engine.log_likelihood_ratio(t, n, theta)


# Internal name used by the calibration engine.
log_mixture = ev_log_mixture

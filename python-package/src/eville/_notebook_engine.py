"""Stable full-line invariant t-mixtures for the publication simulations.

Gaussian priors are integrated analytically in effect size. Other priors
are genuine normal-scale mixtures, integrated on the log-scale axis with
quadrature refinement and explicit bounds for both omitted tails.
Refinement checks are numerical diagnostics, not interval arithmetic.
"""
from dataclasses import dataclass
from functools import lru_cache
import math
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import exp1, gammaln, log_ndtr, logsumexp, roots_legendre

LOG2PI = math.log(2.0 * math.pi)
MIXTURE_RTOL = 2e-9


@dataclass(frozen=True)
class Prior:
    family: str
    location: float = 0.0
    scale: float = 1.0

    def __post_init__(self):
        if self.family not in ('Gaussian', 'Cauchy', 'Student-t5', 'Laplace', 'Horseshoe'):
            raise ValueError('Unknown prior family: ' + str(self.family))
        if not math.isfinite(self.location) or not math.isfinite(self.scale) or self.scale <= 0:
            raise ValueError('The prior requires finite location and positive finite scale.')


def _check_n(n):
    if isinstance(n, (bool, np.bool_)) or not np.isscalar(n) or not np.isfinite(n) or int(n) != n or n < 2:
        raise ValueError('n must be an integer at least 2.')
    return int(n)


def _invariant(t, n):
    n = _check_n(n)
    t = float(t)
    if not math.isfinite(t):
        raise ValueError('Supply a finite t statistic with positive sample variance.')
    h = math.hypot(math.sqrt(n - 1), t)
    u = math.sqrt(n) * (t / h)
    log_d = math.log(n) + math.log(n - 1) - 2 * math.log(h)
    return u, log_d


def _log_chi_constant(n):
    return (n / 2 - 1) * math.log(2) + gammaln(n / 2)


def _log_j_positive(z, n):
    """log integral r^(n-1) exp(-(r-z)^2/2) dr, z >= 0."""
    z = np.asarray(z, dtype=float)
    log_z = np.full_like(z, -np.inf)
    np.log(z, out=log_z, where=z > 0)
    previous = 0.5 * LOG2PI + log_ndtr(z)
    current = np.logaddexp(-0.5 * z * z, log_z + previous)
    for k in range(1, n - 1):
        previous, current = current, np.logaddexp(log_z + current, math.log(k) + previous)
    return current


@lru_cache(maxsize=512)
def _gamma_rule(n, order):
    """Normalized Gauss-Laguerre rule for Q ~ Gamma(n,1)."""
    j = np.arange(order, dtype=float)
    diagonal = 2 * j + n
    k = np.arange(1, order, dtype=float)
    off_diagonal = np.sqrt(k * (k + n - 1))
    nodes, vectors = eigh_tridiagonal(diagonal, off_diagonal)
    weights = vectors[0] ** 2
    good = weights > 0
    return nodes[good], np.log(weights[good])


def _log_i2_negative(a):
    """log integral r exp(-r^2/2-a*r) dr, avoiding cancellation."""
    a = np.asarray(a, dtype=float)
    result = np.empty_like(a)
    small = a < 10
    if np.any(small):
        aa = a[small]
        log_a = np.full_like(aa, -np.inf)
        np.log(aa, out=log_a, where=aa > 0)
        log_product = log_a + 0.5 * LOG2PI + 0.5 * aa * aa + log_ndtr(-aa)
        result[small] = np.log(-np.expm1(log_product))
    if np.any(~small):
        aa = a[~small]
        inv = 1 / (aa * aa)
        term = np.ones_like(aa)
        total = term.copy()
        for j in range(1, 60):
            term *= -(2 * j + 1) * inv
            total += term
            if np.all(np.abs(term) <= 2e-15 * np.abs(total)):
                break
        if np.any(total <= 0) or np.any(np.abs(term) > 2e-13 * np.abs(total)):
            raise FloatingPointError('The negative n=2 radial expansion failed its remainder check.')
        result[~small] = -2 * np.log(aa) + np.log(total)
    return result


def _log_i_negative(a, n, order=128):
    """log integral r^(n-1) exp(-r^2/2-a*r) dr, a >= 0."""
    a = np.asarray(a, dtype=float)
    if n == 2:
        return _log_i2_negative(a)
    nodes, log_weights = _gamma_rule(n, order)
    shape = a.shape
    aa = a.reshape(-1)
    out = np.empty_like(aa)
    # Chunking bounds temporary memory when refining many local scales.
    for start in range(0, aa.size, 2048):
        part = aa[start:start + 2048]
        rate = np.maximum(part, math.sqrt(n))
        values = ((1 - part[:, None] / rate[:, None]) * nodes[None, :]
                  - 0.5 * (nodes[None, :] / rate[:, None]) ** 2
                  + log_weights[None, :])
        out[start:start + part.size] = (gammaln(n) - n * np.log(rate)
                                               + logsumexp(values, axis=1))
    return out.reshape(shape)


def log_A_n(z, n):
    """Stable log E exp(z R), R ~ chi_n, including large positive z."""
    n = _check_n(n)
    z = np.asarray(z, dtype=float)
    if np.any(~np.isfinite(z)):
        raise ValueError('z must be finite.')
    result = np.empty_like(z)
    positive = z >= 0
    result[positive] = (_log_j_positive(z[positive], n)
                        + 0.5 * z[positive] ** 2 - _log_chi_constant(n))
    if np.any(~positive):
        coarse = _log_i_negative(-z[~positive], n, 64)
        fine = _log_i_negative(-z[~positive], n, 128)
        if np.max(np.abs(coarse - fine)) > 2e-10:
            refined = _log_i_negative(-z[~positive], n, 256)
            if np.max(np.abs(refined - fine)) > 2e-10:
                raise FloatingPointError('Negative radial quadrature failed order refinement.')
            fine = refined
        result[~positive] = fine - _log_chi_constant(n)
    return float(result) if result.ndim == 0 else result


def log_likelihood_ratio(t, n, theta):
    """Full signed-effect invariant likelihood ratio, computed in log space."""
    n = _check_n(n)
    u, log_d = _invariant(t, n)
    theta = np.asarray(theta, dtype=float)
    if np.any(~np.isfinite(theta)):
        raise ValueError('Effect sizes must be finite.')
    z = u * theta
    out = np.empty_like(theta)
    positive = z >= 0
    out[positive] = (-0.5 * np.exp(log_d) * theta[positive] ** 2
                     + _log_j_positive(z[positive], n) - _log_chi_constant(n))
    if np.any(~positive):
        a = -z[~positive]
        coarse = _log_i_negative(a, n, 64)
        fine = _log_i_negative(a, n, 128)
        if np.max(np.abs(coarse - fine)) > 2e-10:
            refined = _log_i_negative(a, n, 256)
            if np.max(np.abs(refined - fine)) > 2e-10:
                raise FloatingPointError('Radial likelihood quadrature did not converge.')
            fine = refined
        out[~positive] = -0.5 * n * theta[~positive] ** 2 + fine - _log_chi_constant(n)
    out[theta == 0] = 0
    if np.any(~np.isfinite(out)):
        raise FloatingPointError('Likelihood ratio is outside numerical range.')
    return float(out) if out.ndim == 0 else out


def log_lebesgue(t, n):
    n = _check_n(n)
    _, log_d = _invariant(t, n)
    return 0.5 * math.log(2 * math.pi / n) + 0.5 * n * (math.log(n) - log_d)


def _log_gaussian_from_log_sd(t, n, location, log_sd, radial_order=128):
    u, log_d = _invariant(t, n)
    log_sd = np.asarray(log_sd, dtype=float)
    log_a = np.logaddexp(0, math.log(n) + 2 * log_sd)
    log_D = np.logaddexp(0, log_d + 2 * log_sd)
    z = location * u * np.exp(-0.5 * (log_a + log_D))
    prefactor = 0.5 * (n - 1) * log_a - 0.5 * n * log_D
    if location * u >= 0:
        # Combine the large opposing exponent terms algebraically, before
        # evaluation: -n*m^2/(2a)+z^2/2 = -m^2*d/[2(1+s^2*d)].
        result = (prefactor - 0.5 * location * location * np.exp(log_d - log_D)
                  + _log_j_positive(z, n) - _log_chi_constant(n))
    else:
        result = (prefactor - 0.5 * n * location * location * np.exp(-log_a)
                  + _log_i_negative(-z, n, radial_order) - _log_chi_constant(n))
    return result


def log_gaussian_mixture(t, n, location=0.0, scale=1.0):
    n = _check_n(n)
    Prior('Gaussian', float(location), float(scale))
    result = _log_gaussian_from_log_sd(t, n, float(location), math.log(scale), 128)
    if location * float(t) < 0 and n > 2:
        refined = _log_gaussian_from_log_sd(t, n, float(location), math.log(scale), 256)
        if abs(float(refined - result)) > 2e-10:
            raise FloatingPointError('Gaussian mixture radial quadrature failed refinement.')
        result = refined
    if not np.isfinite(result):
        raise FloatingPointError('Gaussian mixture is outside numerical range.')
    return float(result)


def _log_scale_weight(x, family):
    """Density of X=log(lambda), where theta|lambda~N(m,s^2 lambda^2)."""
    x = np.asarray(x, dtype=float)
    with np.errstate(over='ignore'):
        if family == 'Cauchy':
            return 0.5 * math.log(2 / math.pi) - x - 0.5 * np.exp(-2 * x)
        if family == 'Student-t5':
            a = 2.5
            return math.log(2) + a * math.log(a) - gammaln(a) - 5 * x - a * np.exp(-2 * x)
        if family == 'Laplace':
            return 2 * x - 0.5 * np.exp(2 * x)
        if family == 'Horseshoe':
            return math.log(2 / math.pi) + x - np.logaddexp(0, 2 * x)
    raise ValueError('No scale mixture for this family.')


def horseshoe_logpdf(x, location=0.0, scale=1.0):
    """Genuine N(m,s^2 lambda^2), lambda~half-Cauchy(0,1), density."""
    Prior('Horseshoe', float(location), float(scale))
    z = (np.asarray(x, dtype=float) - location) / scale
    a = 0.5 * z * z
    result = np.empty_like(a)
    zero = a == 0
    regular = (a < 50) & ~zero
    result[zero] = np.inf
    result[regular] = a[regular] + np.log(exp1(a[regular]))
    large = a >= 50
    if np.any(large):
        aa = a[large]
        term = np.ones_like(aa)
        total = term.copy()
        for k in range(1, 25):
            term *= -k / aa
            total += term
            if np.all(np.abs(term) < 2e-15 * np.abs(total)):
                break
        result[large] = -np.log(aa) + np.log(total)
    result -= math.log(scale) + math.log(math.pi) + 0.5 * LOG2PI
    return float(result) if result.ndim == 0 else result


@lru_cache(maxsize=8)
def _legendre(order):
    nodes, weights = roots_legendre(order)
    return nodes, np.log(weights)


def _scale_integral(t, n, prior, lower, upper, order):
    # Every interval is short enough to resolve both broad half-Cauchy scale
    # tails and the narrow variance-mixture saddle for shifted Laplace priors.
    count = max(1, math.ceil((upper - lower) / 0.5))
    edges = np.linspace(lower, upper, count + 1)
    half = np.diff(edges) / 2
    mid = (edges[:-1] + edges[1:]) / 2
    q, log_w = _legendre(order)
    x = (mid[:, None] + half[:, None] * q).ravel()
    weights = (np.log(half[:, None]) + np.broadcast_to(log_w, (count, order))).ravel()
    log_density = _log_scale_weight(x, prior.family)
    live = np.isfinite(log_density)
    value = _log_gaussian_from_log_sd(t, n, prior.location,
                                      math.log(prior.scale) + x[live],
                                      radial_order=max(64, 4 * order))
    answer = logsumexp(value + log_density[live] + weights[live])
    if not np.isfinite(answer):
        raise FloatingPointError('The normal-scale-mixture integral is not finite.')
    return float(answer), int(np.sum(live))


def _scale_tail_logs(t, n, prior, lower, upper):
    """Upper bounds on omitted lower/upper log-scale mixture contributions."""
    _, log_d = _invariant(t, n)
    # L(theta) <= C exp(-k theta^2), hence L(theta) <= C.
    log_C = 0.5 * n * (math.log(2 * n) - log_d)
    family = prior.family
    if family == 'Horseshoe':
        lower_mass = min(0.0, math.log(2 / math.pi) + lower)
    elif family == 'Laplace':
        lower_mass = min(0.0, 2 * lower - math.log(2))
    elif family == 'Cauchy':
        if -lower > 350:
            lower_mass = -math.inf
        else:
            lower_mass = min(0.0, math.log(2) + float(log_ndtr(-math.exp(-lower))))
    elif family == 'Student-t5':
        # Chernoff bound for Gamma(a,1) exceeding a exp(-2 lower).
        if lower >= 0:
            lower_mass = 0.0
        elif -2 * lower > 700:
            lower_mass = -math.inf
        else:
            a = 2.5
            lower_mass = a * (1 - 2 * lower - math.exp(-2 * lower))
    else:
        raise ValueError('Unknown scale-mixture family.')
    # phi(theta;m,sd) <= 1/(sqrt(2pi)*sd), so integrate the Lebesgue
    # likelihood multiplied by an analytic scale-tail moment bound.
    log_factor = log_lebesgue(t, n) - 0.5 * LOG2PI - math.log(prior.scale)
    if family == 'Horseshoe':
        upper_moment = -math.log(math.pi) - 2 * upper
    elif family == 'Cauchy':
        upper_moment = 0.5 * math.log(2 / math.pi) - math.log(2) - 2 * upper
    elif family == 'Student-t5':
        a = 2.5
        upper_moment = math.log(2) + a * math.log(a) - gammaln(a) - math.log(6) - 6 * upper
    else:
        upper_moment = (-math.inf if upper > 350 else
                        0.5 * LOG2PI + float(log_ndtr(-math.exp(upper))))
    return log_C + lower_mass, log_factor + upper_moment


@dataclass(frozen=True)
class MixtureDiagnostic:
    log_value: float
    relative_refinement_change: float
    relative_omitted_tail_bound: float
    log_scale_lower: float
    log_scale_upper: float
    order: int
    nodes: int


@lru_cache(maxsize=50000)
def log_mixture_with_diagnostics(t, n, prior, rtol=MIXTURE_RTOL):
    n = _check_n(n)
    t = float(t)
    _invariant(t, n)
    if not isinstance(prior, Prior) or not 0 < rtol < 0.01:
        raise ValueError('Supply a Prior and a relative tolerance in (0,0.01).')
    if prior.family == 'Gaussian':
        return MixtureDiagnostic(log_gaussian_mixture(t, n, prior.location, prior.scale),
                                  0.0, 0.0, math.nan, math.nan, 0, 0)
    lower, upper = -6.0, 6.0
    for attempt in range(12):
        coarse, _ = _scale_integral(t, n, prior, lower, upper, 16)
        fine, nodes = _scale_integral(t, n, prior, lower, upper, 32)
        discrepancy = abs(math.expm1(coarse - fine))
        order = 32
        if not math.isfinite(discrepancy) or discrepancy > rtol / 4:
            refined, nodes = _scale_integral(t, n, prior, lower, upper, 64)
            discrepancy = abs(math.expm1(fine - refined))
            fine, order = refined, 64
            if not math.isfinite(discrepancy) or discrepancy > rtol / 4:
                raise FloatingPointError('Normal-scale-mixture quadrature failed order refinement.')
        log_low, log_high = _scale_tail_logs(t, n, prior, lower, upper)
        log_tail = float(np.logaddexp(log_low, log_high))
        if log_tail - fine <= math.log(rtol / 4):
            return MixtureDiagnostic(fine, discrepancy, math.exp(log_tail - fine),
                                      lower, upper, order, nodes)
        target = fine + math.log(rtol / 16)
        _, log_d = _invariant(t, n)
        log_C = 0.5 * n * (math.log(2 * n) - log_d)
        if log_low > target:
            if prior.family == 'Horseshoe':
                lower = min(lower - 2, target - log_C - math.log(2 / math.pi))
            elif prior.family == 'Laplace':
                lower = min(lower - 2, 0.5 * (target - log_C + math.log(2)))
            else:
                lower -= 2
        if log_high > target:
            slope = 6 if prior.family == 'Student-t5' else 2
            upper += max(2, (log_high - target) / slope)
        if not math.isfinite(lower + upper) or upper - lower > 4000:
            raise FloatingPointError('Scale separation exceeds the supported numerical range.')
    raise FloatingPointError('Normal-scale-mixture omitted-tail bounds did not converge.')


def log_mixture(t, n, prior, rtol=MIXTURE_RTOL):
    """Scalar log evidence for a fixed full-line location-scale prior."""
    return log_mixture_with_diagnostics(float(t), _check_n(n), prior, rtol).log_value
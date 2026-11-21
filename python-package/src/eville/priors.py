"""Fixed effect-size priors, matching the R package and simulation notebook."""
from dataclasses import dataclass, field
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import betaln, exp1, hyp2f1, log_ndtr
from scipy import stats

_FAMILIES = ('gaussian', 'cauchy', 'student_t5', 'laplace', 'horseshoe', 'lebesgue')
_ALIASES = {'student_t': 'student_t5', 'student-t5': 'student_t5', 't5': 'student_t5'}


def _finite_number(value, name):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(f'{name} must be a finite real number.')
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError):
        raise ValueError(f'{name} must be a finite real number.') from None
    if not math.isfinite(value):
        raise ValueError(f'{name} must be a finite real number.')
    return value


def _log_atan_exp(x):
    if x < -20:
        return x
    if x < 0:
        return math.log(math.atan(math.exp(x)))
    return math.log(math.atan2(1.0, math.exp(-x)))


def _horseshoe_log_sf(log_y):
    """Stable positive tail of Z*lambda; integrate on the half-normal scale."""
    if log_y < -50:
        # The omitted central mass is smaller than double precision here.
        return -math.log(2)
    if log_y < -5:
        # Near zero, the survival integral has a boundary layer of width y.
        # Integrate the small central mass on x=y*r**4 instead, resolving
        # the logarithmic density spike without subtracting two CDFs.
        y = math.exp(log_y)
        anchor = float(_horseshoe_logpdf(y, 0.0, 1.0))
        def central_integrand(r):
            if r == 0:
                return 0.0
            return math.exp(float(_horseshoe_logpdf(y*r**4, 0.0, 1.0))-anchor)*4*r**3
        value, error = quad(central_integrand, 0, 1, epsrel=1e-11, epsabs=0, limit=200)
        if value <= 0 or error > 1e-9*value:
            raise FloatingPointError('Horseshoe central mass did not converge.')
        central_mass = math.exp(log_y+anchor+math.log(value))
        return -math.log(2)+math.log1p(-2*central_mass)
    anchor = _log_atan_exp(-log_y)
    def integrand(z):
        if z == 0:
            return 0.0
        return math.exp(-0.5*z*z + _log_atan_exp(math.log(z)-log_y) - anchor)
    value, error = quad(integrand, 0, np.inf, epsrel=1e-11, epsabs=0, limit=1000)
    if not math.isfinite(value) or value <= 0 or error > 1e-9 * value:
        raise FloatingPointError('Horseshoe half-line normalization did not converge.')
    return min(-math.log(2), 0.5*math.log(2/math.pi)-math.log(math.pi)+anchor+math.log(value))


def _student_log_sf(y):
    if y < math.sqrt(5):
        return float(stats.t.logsf(y, 5))
    # The incomplete-beta representation evaluated before the tail underflows.
    log_z = math.log(5) - float(np.logaddexp(math.log(5), 2*math.log(y)))
    z = math.exp(log_z)
    return (-math.log(2) + 2.5*log_z - math.log(2.5)
            - betaln(2.5, 0.5) + math.log(float(hyp2f1(2.5, 0.5, 3.5, z))))


def _log_masses(family, location, scale):
    if location == 0:
        return -math.log(2), -math.log(2)
    y = abs(location)/scale
    if family == 'gaussian':
        small = float(log_ndtr(-y))
    elif family == 'cauchy':
        small = _log_atan_exp(math.log(scale)-math.log(abs(location))) - math.log(math.pi)
    elif family == 'student_t5':
        small = _student_log_sf(y)
    elif family == 'laplace':
        small = -math.log(2)-y
    else:
        small = _horseshoe_log_sf(math.log(abs(location))-math.log(scale))
    if not math.isfinite(small):
        raise ValueError('A prior half-line mass is outside numerical range.')
    large = math.log(-math.expm1(small))
    return (large, small) if location > 0 else (small, large)


@dataclass(frozen=True)
class Prior:
    """A fixed prior on the standardized effect (mu-mu0)/sigma.

    ``scale`` defaults to 1/sqrt(2) for Cauchy/Rouder and 1 for the
    other proper families. One-sided distributions are conditioned on
    the chosen half-line. Location and scale describe the distribution
    *before* conditioning. Lebesgue measure has intensity one on its support.
    """
    family: str = 'gaussian'
    location: float = 0.0
    scale: float | None = None
    alternative: str = 'two.sided'
    proper: bool = field(init=False)
    log_positive_mass: float = field(init=False, repr=False)
    log_negative_mass: float = field(init=False, repr=False)

    def __post_init__(self):
        if not isinstance(self.family, str):
            raise ValueError('family must be a prior-family name.')
        family = _ALIASES.get(self.family.lower(), self.family.lower())
        if family not in _FAMILIES + ('rouder',):
            raise ValueError(f'Unknown prior family {self.family!r}.')
        if self.alternative not in ('two.sided', 'less', 'greater'):
            raise ValueError('alternative must be two.sided, less, or greater.')
        location = _finite_number(self.location, 'location')
        if family == 'lebesgue':
            if location != 0 or self.scale is not None:
                raise ValueError('Omit location and scale for Lebesgue measure.')
            scale = None
            positive = negative = math.nan
        else:
            if family == 'rouder':
                if location != 0:
                    raise ValueError('For a shifted Cauchy use family="cauchy".')
                family = 'cauchy'
            scale = (1/math.sqrt(2) if family == 'cauchy' else 1.0) if self.scale is None else _finite_number(self.scale, 'scale')
            if scale <= 0:
                raise ValueError('scale must be positive.')
            if location != 0 and scale <= 16*np.finfo(float).eps*abs(location):
                raise ValueError('The prior scale is below numerical resolution at this location.')
            if not math.isfinite(location/scale):
                raise ValueError('The location/scale ratio exceeds numerical range.')
            positive, negative = _log_masses(family, location, scale)
        object.__setattr__(self, 'family', family)
        object.__setattr__(self, 'location', location)
        object.__setattr__(self, 'scale', scale)
        object.__setattr__(self, 'proper', family != 'lebesgue')
        object.__setattr__(self, 'log_positive_mass', positive)
        object.__setattr__(self, 'log_negative_mass', negative)


def ev_prior(family='gaussian', location=0.0, scale=None, alternative='two.sided'):
    """Construct a :class:`Prior`; fix it independently of the target data."""
    return Prior(family, location, scale, alternative)


def _log_abs_difference(theta, location):
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        difference = np.abs(theta-location)
        log_difference = np.log(difference)
        overflow = np.isinf(difference) & np.isfinite(theta)
        return np.where(overflow, np.logaddexp(np.log(np.abs(theta)), math.log(abs(location)) if location else -np.inf), log_difference)


def _horseshoe_logpdf(theta, location, scale):
    theta = np.asarray(theta, dtype=float)
    log_a = 2*(_log_abs_difference(theta, location)-math.log(scale))-math.log(2)
    out = np.empty_like(log_a)
    zero = np.isneginf(log_a)
    tiny = np.isfinite(log_a) & (log_a < -30)
    regular = (log_a >= -30) & (log_a < math.log(50))
    large = log_a >= math.log(50)
    out[zero] = np.inf
    out[tiny] = np.log(-0.5772156649015329-log_a[tiny])
    a = np.exp(log_a[regular])
    out[regular] = a + np.log(exp1(a))
    inv = np.exp(-log_a[large])
    term = np.ones_like(inv)
    total = term.copy()
    for k in range(1, 26):
        term *= -k*inv
        total += term
        if np.all(np.abs(term) <= 2e-15*np.abs(total)):
            break
    out[large] = -log_a[large]+np.log(total)
    return out-math.log(scale)-math.log(math.pi)-0.5*math.log(2*math.pi)


def _logpdf(theta, family, location, scale):
    theta = np.asarray(theta, dtype=float)
    if family == 'horseshoe':
        return _horseshoe_logpdf(theta, location, scale)
    with np.errstate(over='ignore'):
        if family == 'gaussian':
            return stats.norm.logpdf(theta, loc=location, scale=scale)
        if family == 'cauchy':
            log_z = _log_abs_difference(theta, location)-math.log(scale)
            return -math.log(math.pi)-math.log(scale)-np.logaddexp(0, 2*log_z)
        if family == 'student_t5':
            log_z = _log_abs_difference(theta, location)-math.log(scale)
            constant = math.lgamma(3)-math.lgamma(2.5)-0.5*math.log(5*math.pi)
            return constant-math.log(scale)-3*np.logaddexp(0, 2*log_z-math.log(5))
        if family == 'laplace':
            return -math.log(2)-math.log(scale)-np.abs((theta-location)/scale)
    raise ValueError('Unknown proper prior.')


def ev_prior_density(theta, prior=None, log=False):
    """Evaluate a normalized prior density (Lebesgue: intensity, not density)."""
    prior = ev_prior() if prior is None else prior
    if not isinstance(prior, Prior):
        raise TypeError('prior must be constructed with ev_prior().')
    if not isinstance(log, (bool, np.bool_)):
        raise ValueError('log must be a boolean.')
    values = np.asarray(theta, dtype=float)
    if np.any(~np.isfinite(values)):
        raise ValueError('theta must contain finite real values.')
    inside = np.ones(values.shape, dtype=bool)
    mass = 0.0
    if prior.alternative == 'greater':
        inside = values > 0
        mass = prior.log_positive_mass
    elif prior.alternative == 'less':
        inside = values < 0
        mass = prior.log_negative_mass
    answer = np.full(values.shape, -np.inf)
    if not prior.proper:
        answer[inside] = 0.0
    elif np.any(inside):
        answer[inside] = _logpdf(values[inside], prior.family, prior.location, prior.scale)-mass
    if not log:
        with np.errstate(over='ignore'):
            answer = np.exp(answer)
    return float(answer) if answer.ndim == 0 else answer

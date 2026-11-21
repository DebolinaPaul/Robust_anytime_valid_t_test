"""One-sided invariant Gaussian t-mixtures and extended-Ville calibration.

Use X_i=-Chronos_i, so the scientific alternative is delta=E[X]/sd(X)>0.
The null is E[X]<=0, under the independent Gaussian working model.  Priors
are proper location-scale densities fitted on SOURCE effects alone, then
conditioned on delta>0, or the separately implemented improper flat measure.
The latter uses only extended Ville. n counts target cell-line observations.

For n>=2 and t=sqrt(n)*mean(X)/sd(X),
 L_n(delta)=nct.pdf(t,n-1,sqrt(n)*delta)/t.pdf(t,n-1).
This is the likelihood ratio on the nested positive-scale-invariant
filtration.  For delta>0 it is increasing in t, so the boundary null mean=0
is least favourable.  Coupling X_i=mu+sigma*Z_i shows joint rejection events
are increasing in mu (sample SDs do not change with mu).  Consequently the
boundary calibration also controls the composite null mu<=0.

The exact theory assumes an independently fixed prior, target lineage and
sampling order.  Searching target results is exploratory and is not made
confirmatory by this calibration.  Chronos scores are observational derived
measurements; Gaussianity/independence remain working assumptions.

The fit_prior function supplies median/IQR initial values for optimization
and the explicitly labelled horseshoe fit. The fitting module uses SciPy MLE
for the four other full families. The exact normal--half-Cauchy horseshoe
marginal has no ordinary finite free-location MLE.

Numerical quadrature is approximate; compare quadrature orders and independent
integrals before interpreting values near alpha. Dependencies: numpy, scipy.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import exp1, erfcx, gammaincc, gammaln, log_ndtr, logsumexp, ndtr, owens_t
from scipy.stats import norm, cauchy, laplace, nct, t as student_t
from improper_flat import FlatPrior, log_flat_mixture, flat_calibration, flat_log_threshold

FAMILIES = ('Gaussian', 'Cauchy', 'Laplace', 'Student-t5', 'Horseshoe')
_LOG2PI = math.log(2*math.pi)


def _family(name):
    name = {'Normal': 'Gaussian', 'Student-t_5': 'Student-t5', 't5': 'Student-t5',
            'Cauchy (Rouder)': 'Cauchy', 'Rouder': 'Cauchy'}.get(name, name)
    if name not in FAMILIES:
        raise ValueError(f'Unknown family {name}; choose from {FAMILIES}')
    return name


def _log_scaled_e1(a):
    """log(exp(a)*E1(a)), avoiding overflow and logarithmic singularity."""
    a = np.maximum(np.asarray(a, float), np.finfo(float).tiny)
    out = np.empty_like(a)
    small = a < 50.
    out[small] = a[small]+np.log(exp1(a[small]))
    if np.any(~small):
        b = a[~small]
        term = np.ones_like(b)
        value = term.copy()
        for k in range(1, 18):
            term *= -k/b
            value += term
        out[~small] = -np.log(b)+np.log(value)
    return out


def _hs_logpdf(z):
    z = np.asarray(z, float)
    value = _log_scaled_e1(.5*z**2)-math.log(math.pi)-.5*_LOG2PI
    # The exact horseshoe density diverges at its location. Quadrature splits
    # there and uses interior nodes, so this measure-zero value changes no
    # mixture integral; it must not be a finite artificial likelihood peak.
    return np.where(z == 0, np.inf, value)


@lru_cache(maxsize=2048)
def _hs_sf_positive(z):
    """Horseshoe survival probability for finite z >= 0."""
    z = float(z)
    if not math.isfinite(z) or z < 0:
        raise ValueError("z must be finite and nonnegative")
    if z == 0:
        return 0.5

    if z < 1:
        def integrand(angle): return ndtr(-z / math.tan(angle))
        upper = math.pi / 2
    else:
        # Exact substitution u = z/lambda avoids a narrow endpoint.
        def integrand(u): return ndtr(-u) / (1 + (u / z)**2)
        upper = np.inf

    result = quad(
        integrand, 0.0, upper,
        epsabs=2e-13, epsrel=2e-12,
        limit=250, full_output=1
    )
    value, error = result[:2]
    tolerance = max(2e-13, 2e-12 * abs(value))

    if (len(result) != 3
            or not np.isfinite([value, error]).all()
            or error > tolerance):
        raise FloatingPointError(
            "Horseshoe survival integration did not converge"
        )

    # Divide before multiplying to avoid overflow in pi*z.
    answer = (2 / math.pi) * (value if z < 1 else value / z)
    if not 0 < answer <= 0.5:
        raise FloatingPointError(
            "Horseshoe survival probability is out of range"
        )
    return float(answer)


@lru_cache(maxsize=16)
def _q75(family):
    if family == 'Horseshoe':
        return brentq(lambda x: _hs_sf_positive(float(x))-.25, 1e-8, 20., xtol=1e-13)
    return {'Gaussian': norm, 'Cauchy': cauchy, 'Laplace': laplace,
            'Student-t5': student_t(5)}[family].ppf(.75)


@dataclass(frozen=True)
class Prior:
    family: str
    location: float
    scale: float

    def __post_init__(self):
        object.__setattr__(self, 'family', _family(self.family))
        if not np.isfinite(self.location) or not np.isfinite(self.scale) or self.scale <= 0:
            raise ValueError(
                'Finite location and positive finite scale required')
        z0 = -self.location/self.scale
        if self.family == 'Horseshoe':
            sf = _hs_sf_positive(float(abs(z0)))
            log_z = math.log(sf if z0 >= 0 else 1-sf)
        else:
            dist = {'Gaussian': norm, 'Cauchy': cauchy, 'Laplace': laplace,
                    'Student-t5': student_t(5)}[self.family]
            log_z = float(dist.logsf(z0))
        if not np.isfinite(log_z):
            raise ValueError(
                'Positive-half prior mass is below floating-point resolution')
        object.__setattr__(self, 'log_positive_mass', log_z)

    @property
    def positive_mass(self):
        return math.exp(self.log_positive_mass)

    def untruncated_logpdf(self, delta):
        z = (np.asarray(delta, float)-self.location)/self.scale
        if self.family == 'Gaussian':
            ans = -.5*_LOG2PI-.5*z*z
        elif self.family == 'Cauchy':
            ans = -math.log(math.pi)-np.log1p(z*z)
        elif self.family == 'Laplace':
            ans = -math.log(2)-np.abs(z)
        elif self.family == 'Student-t5':
            ans = gammaln(3)-gammaln(2.5)-.5 * \
                math.log(5*math.pi)-3*np.log1p(z*z/5)
        else:
            ans = _hs_logpdf(z)
        return ans-math.log(self.scale)

    def logpdf(self, delta):
        delta = np.asarray(delta, float)
        return np.where(delta > 0, self.untruncated_logpdf(delta)-self.log_positive_mass, -np.inf)

    def pdf(self, delta):
        return np.exp(self.logpdf(delta))


def fit_prior(effects, family, min_scale=1e-4):
    """Median/IQR fit of full family, followed by positive-half truncation.

    'effects' contains one estimate mean(-Chronos)/sd(Chronos) per SOURCE
    lineage.  Each lineage has equal weight.  Quantile matching avoids the
    unbounded location likelihood of an exact horseshoe density.
    """
    x = np.asarray(effects, float)
    if x.ndim != 1 or len(x) < 4 or not np.isfinite(x).all():
        raise ValueError(
            'At least four finite source-lineage effects are required')
    family = _family(family)
    q25, q75 = np.quantile(x, [.25, .75])
    scale = max(float((q75-q25)/(2*_q75(family))), float(min_scale))
    return Prior(family, float(np.median(x)), scale)


@lru_cache(maxsize=32)
def _legendre(order):
    return np.polynomial.legendre.leggauss(order)


def _log_chi_constant(n):
    return (.5*n-1)*math.log(2)+gammaln(.5*n)


def _log_j_positive(z, n):
    """log integral r^(n-1)*exp(-(r-z)^2/2) dr for z>=0.

    Integer-n positive recurrence is stable, including in extreme tails.
    """
    z = np.asarray(z, float)
    j1 = .5*_LOG2PI+log_ndtr(z)
    log_z = np.full_like(z, -np.inf)
    np.log(z, out=log_z, where=z > 0)
    j2 = np.logaddexp(-.5*z*z, log_z+j1)
    if n == 2:
        return j2
    for k in range(1, n-1):
        j1, j2 = j2, np.logaddexp(log_z+j2, math.log(k)+j1)
    return j2


def _log_j_negative(z, n):
    """Positive modal quadrature for z<=0; no subtractive recurrences."""
    z = np.asarray(z, float)
    if n == 2:
        a = -z
        factor = np.empty_like(a)
        small = a < 30.
        factor[small] = 1-a[small] * \
            math.sqrt(math.pi/2)*erfcx(a[small]/math.sqrt(2))
        if np.any(~small):
            inv = 1/(a[~small]*a[~small])
            term = inv.copy()
            value = term.copy()
            for k in range(1, 15):
                term *= -(2*k+1)*inv
                value += term
            factor[~small] = value
        return -.5*a*a+np.log(factor)
    flat = z.ravel()
    d = np.hypot(flat, 2*math.sqrt(n-1))
    mode = 2*(n-1)/(d-flat)
    sd = 1/np.sqrt((n-1)/mode**2+1)
    # Mode-centred intervals resolve both a sharp mode and the r=0 boundary.
    offsets = np.array([-40, -24, -16, -12, -8, -6, -4, -
                       2, -1, 0, 1, 2, 4, 6, 8, 12, 16, 24, 40.])
    breaks = np.maximum(0., mode[:, None]+sd[:, None]*offsets)
    half = np.diff(breaks, axis=1)/2
    mid = (breaks[:, 1:]+breaks[:, :-1])/2
    gx, gw = _legendre(16)
    r = mid[:, :, None]+half[:, :, None]*gx
    with np.errstate(divide='ignore', invalid='ignore'):
        val = ((n-1)*np.log(r)-.5*(r-flat[:, None, None])**2
               + np.log(half[:, :, None])+np.log(gw))
    val = np.where(half[:, :, None] > 0, val, -np.inf)
    return logsumexp(val, axis=(1, 2)).reshape(z.shape)


def log_likelihood_ratio(t, n, delta):
    if (isinstance(n, (bool, np.bool_)) or not np.isscalar(n)
            or not np.isfinite(n) or int(n) != n or n < 2):
        raise ValueError('n must be a finite integer >= 2')
    n = int(n)
    if not np.isscalar(t) or not np.isfinite(t):
        raise ValueError('t must be a finite scalar')
    delta = np.asarray(delta, float)
    if not np.isfinite(delta).all() or np.any(delta < 0):
        raise ValueError('This one-sided engine requires finite delta >= 0')
    hyp = math.hypot(math.sqrt(n-1), float(t))
    h = math.sqrt(n)*float(t)/hyp
    v = n*(math.sqrt(n-1)/hyp)**2
    z = h*delta
    log_j = _log_j_positive(z, n) if t >= 0 else _log_j_negative(z, n)
    return -.5*v*delta**2+log_j-_log_chi_constant(n)


def _mixture_nodes(t, n, prior, order):
    hyp = math.hypot(math.sqrt(n-1), float(t))
    sigma_like = hyp/math.sqrt(n*(n-1))
    likelihood_center = max(0., float(t)/math.sqrt(n))
    centers = np.unique([0., max(prior.location, 0.), likelihood_center])
    scales = np.unique([sigma_like, prior.scale])
    offsets = np.r_[-2.**np.arange(6, -4, -1), 0., 2.**np.arange(-3, 7)]
    points = (centers[:, None, None]+scales[None, :, None]
              * offsets[None, None, :]).ravel()
    points = np.r_[0., points[points > 0]]
    if prior.family == 'Horseshoe' and prior.location > 0:
        tiny = prior.scale*10.**np.arange(-12, 0)
        points = np.r_[points, prior.location,
                       prior.location-tiny, prior.location+tiny]
    points = np.unique(points[points >= 0])
    # Bridge enormous scale separations geometrically. This is essential for
    # Cauchy/horseshoe tails when the initial null t statistic is extreme:
    # L_2(delta)*pi(delta) may behave as 1/delta over many decades.
    geometric_low = max(min(scales)*.0625, 1e-12)
    geometric_high = float(points[-1])
    count = max(2, int(math.ceil(2*math.log2(geometric_high/geometric_low)))+1)
    points = np.unique(np.r_[points, np.geomspace(
        geometric_low, geometric_high, count)])
    # Avoid integrating enormous empty prior tails for light families. Keeping
    # the union of likelihood/prior intervals is robust to displaced priors.
    half, mid = np.diff(points)/2, (points[:-1]+points[1:])/2
    gx, gw = _legendre(order)
    delta_grid = mid[:, None]+half[:, None]*gx
    if prior.family == 'Horseshoe' and prior.location > 0:
        # The exact density is infinite at its location. Interior nodes on
        # intervals ending there can round to that endpoint at high orders;
        # move only such nodes one floating-point step into their interval.
        rounded_to_spike = delta_grid == prior.location
        if np.any(rounded_to_spike):
            interior = np.broadcast_to(mid[:, None], delta_grid.shape)
            delta_grid = np.where(rounded_to_spike,
                                  np.nextafter(prior.location, interior), delta_grid)
    delta = delta_grid.ravel()
    log_w = (np.log(half[:, None])+np.log(gw)).ravel()
    return delta, log_w


def log_mixture(t, n, prior, quadrature_order=24):
    if isinstance(prior, FlatPrior):
        return log_flat_mixture(t, n)
    delta, log_w = _mixture_nodes(t, n, prior, int(quadrature_order))
    return float(logsumexp(log_likelihood_ratio(t, n, delta)+prior.logpdf(delta)+log_w))


def _df1_nct_cdf(t_value, nc, order=32):
    """Stable noncentral-t CDF for df=1 and nonnegative noncentrality.

    For t>=0 use the exact positive-term Owen-T identity. For t<0 use
      F(t;1,nc) = 2 integral_0^infinity phi(w) Phi(t*w-nc) dw,
    evaluated with positive quadrature weights. This avoids NaNs from some
    SciPy/Boost noncentral-t CDF implementations in the negative tail.
    """
    nc = np.asarray(nc, float)
    if not np.isfinite(t_value) or not np.isfinite(nc).all() or np.any(nc < 0):
        raise ValueError(
            'Finite t and nonnegative finite noncentralities required')
    if t_value >= 0:
        h = nc/math.hypot(1., t_value)
        answer = ndtr(-h)+2*owens_t(h, t_value)
        # owens_t can return a negative subnormal near underflow. Reject
        # meaningful range violations but tolerate absolute roundoff.
        if (not np.isfinite(answer).all() or np.any(answer < -1e-14)
                or np.any(answer > 1+1e-14)):
            raise FloatingPointError(
                'Owen-T CDF failed beyond roundoff tolerance')
        return np.clip(answer, 0., 1.)
    a = -float(t_value)
    shape = nc.shape
    flat = nc.ravel()
    # Resolve the width at w=0 both for large |t| and large nc. After this
    # rescaling, truncating at x=64 omits a negligible exponentially small
    # tail; independent adaptive integration and order refinement verify it.
    rate = math.hypot(1., a)+a*flat
    edges = np.array([0., .25, .5, 1., 2., 4., 8., 16., 32., 64.])
    half, mid = np.diff(edges)/2, (edges[:-1]+edges[1:])/2
    gx, gw = _legendre(order)
    xx = (mid[:, None]+half[:, None]*gx).ravel()
    log_weights = (np.log(half[:, None])+np.log(gw)).ravel()
    ww = xx[None, :]/rate[:, None]
    log_integrand = (math.log(2)-.5*_LOG2PI-.5*ww*ww
                     + log_ndtr(-flat[:, None]-a*ww)-np.log(rate[:, None])+log_weights[None, :])
    return np.exp(logsumexp(log_integrand, axis=1)).reshape(shape)


def calibration_direct(log_c, prior, n0=2, quadrature_order=48):
    """Independent clipped-expectation check using noncentral-t CDFs.

    If t_star solves M_n0(t_star)=c, monotonicity gives
      g(c)=P_0(T>t_star)+c^{-1} integral P_delta(T<=t_star) pi(delta)d delta.
    This avoids interpolating the initial evidence curve. If the crossing
    lies beyond |t|=sinh(40), return the exact proper-prior upper bound 1/c.
    """
    if n0 != 2:
        raise ValueError(
            'This audited package implements calibration only for n0=2.')
    if isinstance(prior, FlatPrior):
        return flat_calibration(log_c, n0=n0)
    log_c = float(log_c)

    def objective(x):
        return log_mixture(math.sinh(x), n0, prior, quadrature_order)-log_c
    if objective(40.) <= 0:
        return math.exp(min(0., -log_c))
    if objective(-40.) >= 0:
        return 1.
    root = brentq(objective, -40., 40., xtol=1e-10)
    t_star = math.sinh(root)
    delta, log_w = _mixture_nodes(t_star, n0, prior, quadrature_order)
    nc = np.sqrt(n0)*delta
    probs = np.empty_like(nc)
    ordinary = nc <= 50
    if n0 == 2:
        probs = _df1_nct_cdf(t_star, nc)
        ordinary = np.ones_like(nc, dtype=bool)
    else:
        probs[ordinary] = nct.cdf(t_star, n0-1, nc[ordinary])
    if np.any(~ordinary):
        if t_star <= 0:
            probs[~ordinary] = 0.
        else:
            # An independent positive integral avoids extremely slow Boost
            # noncentral-t series at huge noncentrality. For nc>50, the
            # omitted Z<-nc Gaussian tail is below floating-point precision.
            gx, gw = np.polynomial.hermite.hermgauss(32)
            ratio = (nc[~ordinary, None]+np.sqrt(2)*gx)/t_star
            probs[~ordinary] = np.sum(gammaincc((n0-1)/2,
                                                .5*(n0-1)*ratio**2)*gw[None, :], axis=1)/math.sqrt(math.pi)
    if not np.isfinite(probs).all() or np.any(probs < -1e-14) or np.any(probs > 1 + 1e-12):
        raise FloatingPointError(
            'Noncentral-t CDF evaluation failed; cannot calibrate safely')
    probs = np.clip(probs, 0., 1.)
    with np.errstate(divide='ignore'):
        log_cdf_integral = logsumexp(np.log(probs)+prior.logpdf(delta)+log_w)
    p = student_t.sf(t_star, n0-1)+math.exp(log_cdf_integral-log_c)
    # Both exact quantities are bounds. Taking their minimum is justified for
    # the exact integrals, but is NOT a certified remedy for quadrature error.
    return min(1., math.exp(min(0., -log_c)), float(p))


def direct_log_threshold(prior, alpha=.05, n0=2, quadrature_order=48):
    """Numerical proper-prior boundary; analytic flat-prior calibration."""
    if not 0 < alpha < 1 or n0 != 2:
        raise ValueError('Require 0<alpha<1 and n0=2 in this package.')
    if isinstance(prior, FlatPrior):
        return flat_log_threshold(alpha=alpha, n0=n0)
    upper = math.log(1/alpha)
    # When clipping offers no numerically resolvable improvement, the exact
    # proper-prior Ville threshold is valid and avoids a roundoff-only root.
    if calibration_direct(upper, prior, n0, quadrature_order) >= alpha*(1-1e-10):
        return upper
    lower = min(0., upper-1.)
    while calibration_direct(lower, prior, n0, quadrature_order) < alpha:
        lower -= 10
    return float(brentq(lambda x: calibration_direct(x, prior, n0, quadrature_order)-alpha,
                        lower, upper, xtol=2e-10))


def trajectory(x, prior, calibration=None, alpha=.05, n0=2,
               quadrature_order=24, compare_ville=True):
    """Running evidence with stable variance and direct default calibration."""
    if n0 != 2 or not 0 < alpha < 1:
        raise ValueError("Require n0=2 and 0 < alpha < 1")

    if (not np.isfinite(quadrature_order)
            or int(quadrature_order) != quadrature_order
            or quadrature_order < 2):
        raise ValueError("quadrature_order must be an integer >= 2")
    quadrature_order = int(quadrature_order)

    if isinstance(prior, FlatPrior) and compare_ville:
        raise ValueError(
            "Use compare_ville=False for the improper flat prior"
        )

    if calibration is not None:
        if getattr(calibration, "n0", n0) != n0:
            raise ValueError("Calibration and trajectory have different n0")
        if getattr(calibration, "prior", prior) != prior:
            raise ValueError(
                "Calibration and trajectory have different priors")

    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError("Provide at least two finite observations")

    ns = np.arange(2, len(x) + 1)
    ts = np.empty(len(ns))
    origin = float(x[0])
    mean_offset, m2 = 0.0, 0.0

    for n, value in enumerate(x, start=1):
        # Welford applied to offsets improves accuracy for nearly equal data.
        offset = float(value) - origin
        difference = offset - mean_offset
        mean_offset += difference / n
        m2 += difference * (offset - mean_offset)

        if n >= 2:
            variance = m2 / (n - 1)
            if not math.isfinite(variance) or variance <= 0:
                raise ValueError(
                    "Each monitored prefix needs a finite, positive "
                    "variance; check for ties or overflow"
                )
            ts[n - 2] = math.sqrt(n) * (
                (origin + mean_offset) / math.sqrt(variance)
            )

    if not np.isfinite(ts).all():
        raise FloatingPointError(
            "A t statistic exceeded floating-point range"
        )

    log_m = np.array([
        log_mixture(t, int(n), prior, quadrature_order)
        for t, n in zip(ts, ns)
    ])
    if not np.isfinite(log_m).all():
        raise FloatingPointError(
            "Nonfinite mixture evidence; inspect quadrature"
        )
    maxima = np.maximum.accumulate(log_m)

    if calibration is None:
        # Evaluate repeated maxima only once. calibration_direct also
        # dispatches the improper flat prior to its analytic calibration.
        unique, inverse = np.unique(maxima, return_inverse=True)
        values = np.array([
            calibration_direct(
                float(v), prior, n0=2,
                quadrature_order=max(48, quadrature_order)
            )
            for v in unique
        ])
        p_e = values[inverse]
    else:
        p_e = np.asarray(calibration.pvalue(maxima), dtype=float)

    if (p_e.shape != maxima.shape
            or not np.isfinite(p_e).all()
            or np.any((p_e < 0) | (p_e > 1))):
        raise FloatingPointError("Calibration returned invalid p-values")

    p_v = (
        np.exp(np.minimum(0.0, -maxima))
        if compare_ville else np.full_like(maxima, np.nan)
    )

    def first_rejection(p):
        hits = np.flatnonzero(p <= alpha)
        return int(ns[hits[0]]) if len(hits) else None

    return dict(
        n=ns,
        t=ts,
        log_evidence=log_m,
        log_running_evidence=maxima,
        p_eville=p_e,
        p_ville=p_v,
        eville_rejection=first_rejection(p_e),
        ville_rejection=first_rejection(p_v) if compare_ville else None,
        ville_computed=bool(compare_ville)
    )

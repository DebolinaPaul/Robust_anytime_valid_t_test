"""Return source-fitted priors and fixed benchmarks.

    source_effects contains mean(-Chronos)/sd(Chronos) per source lineage.
    source_n contains the corresponding sample sizes in the same order.
    Four families use full-family MLE on the effect estimates.
    Horseshoe uses constrained marginal MLE with noncentral-t sampling.
    Rouder and improper flat are fixed comparisons.
    """
from __future__ import annotations

import math
import warnings

import numpy as np
from scipy import optimize, stats

from engine import FAMILIES, Prior, fit_prior
from improper_flat import FlatPrior


ROUDER_LABEL = 'Rouder half-Cauchy'
IMPROPER_LABEL = 'Improper flat'


def _validated_effects(effects):
    x = np.asarray(effects, dtype=float)
    if x.ndim != 1 or len(x) < 4 or not np.isfinite(x).all():
        raise ValueError(
            'At least four finite source-lineage effects are required')
    if np.ptp(x) <= 0:
        raise ValueError(
            'Source effects must not all coincide: a positive-scale MLE does not exist')
    return x


def _full_nll(prior, x):
    # The exact horseshoe spike gives infinite log likelihood at a data-aligned
    # location, so no finite negative log likelihood should be reported there.
    if prior.family == 'Horseshoe' and np.any(x == prior.location):
        return None
    value = -float(np.sum(prior.untruncated_logpdf(x)))
    return value if np.isfinite(value) else None


def _scipy_mle(x, family, iqr_prior):
    """Fit a full family; multistart numerical fits are checked, not certified global."""
    if family == 'Gaussian':
        location, scale = stats.norm.fit(x)
        details = dict(scipy_function='scipy.stats.norm.fit', optimizer='closed form',
                       optimizer_converged=True, starts_attempted=1, starts_converged=1)
    elif family == 'Laplace':
        location, scale = stats.laplace.fit(x)
        details = dict(scipy_function='scipy.stats.laplace.fit', optimizer='closed form',
                       optimizer_converged=True, starts_attempted=1, starts_converged=1)
    elif family in ('Cauchy', 'Student-t5'):
        # Standardizing the observations improves optimization conditioning;
        # MLE location/scale equivariance maps the answer back exactly.
        center = float(np.mean(x))
        unit = float(np.std(x, ddof=0))
        z = (x-center)/unit
        dist = stats.cauchy if family == 'Cauchy' else stats.t
        fixed = {} if family == 'Cauchy' else {'fdf': 5}
        initial_locations = np.unique(
            np.r_[np.quantile(z, [.25, .5, .75]), 0.])
        initial_scale = iqr_prior.scale/unit
        candidates, failures = [], []

        def scipy_optimizer(func, start, args=(), disp=0):
            # SciPy .fit supplies a loc/scale objective after fixing t df.
            # Optimize log(scale), rather than allow invalid negative scales.
            def objective(parameters):
                if not np.isfinite(parameters).all() or abs(parameters[1]) > 700:
                    return np.inf
                return float(func([parameters[0], math.exp(parameters[1])], *args))

            result = optimize.minimize(
                objective, [start[0], math.log(start[1])], method='Nelder-Mead',
                options={'xatol': 1e-10, 'fatol': 1e-10, 'maxiter': 5000})
            if not result.success or not np.isfinite(result.fun):
                raise RuntimeError(
                    f'SciPy likelihood optimization failed: {result.message}')
            return np.array([result.x[0], math.exp(result.x[1])])

        for initial_location in initial_locations:
            for factor in (.5, 1., 2.):
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter('error', RuntimeWarning)
                        fitted = dist.fit(z, loc=float(initial_location),
                                          scale=float(initial_scale*factor),
                                          method='MLE', optimizer=scipy_optimizer, **fixed)
                    loc_z, scale_z = fitted[-2:]
                    candidate = Prior(family, center+unit*loc_z, unit*scale_z)
                    nll = _full_nll(candidate, x)
                    if nll is None:
                        raise RuntimeError(
                            'Nonfinite likelihood at fitted parameters')
                    candidates.append((nll, candidate))
                except (ValueError, RuntimeError, RuntimeWarning, FloatingPointError) as exc:
                    failures.append(str(exc))
        if not candidates:
            raise RuntimeError(f'Every {family} MLE start failed: {failures}')
        _, chosen = min(candidates, key=lambda pair: pair[0])
        location, scale = chosen.location, chosen.scale
        details = dict(
            scipy_function=('scipy.stats.cauchy.fit' if family == 'Cauchy'
                            else 'scipy.stats.t.fit(fdf=5)'),
            optimizer='multistart Nelder-Mead in location/log-scale through scipy distribution.fit',
            optimizer_converged=True,
            starts_attempted=len(candidates)+len(failures), starts_converged=len(candidates),
            converged_start_nll_spread=float(
                np.ptp([v[0] for v in candidates])),
            optimizer_failures='; '.join(failures),
            global_optimum_certified=False)
    else:
        raise ValueError(
            f'Ordinary location-scale MLE is not implemented for {family}')
    prior = Prior(family, float(location), float(scale))
    nll, baseline_nll = _full_nll(prior, x), _full_nll(iqr_prior, x)
    if nll is None or baseline_nll is None or nll > baseline_nll+1e-7:
        raise RuntimeError(
            f'{family} MLE failed the finite likelihood / IQR baseline check')
    details.update(full_family_nll=nll, iqr_full_family_nll=baseline_nll,
                   nll_improvement_over_iqr=baseline_nll-nll)
    return prior, details


def fit_horseshoe_mml(source_effects, source_n, min_scale=1e-4):
    """Source-only constrained marginal MLE for a full horseshoe family.

    Effects use delta=-theta. Fit on the full real line; Prior subsequently
    conditions the fitted distribution on delta>0 for target testing.
    """
    from scipy.integrate import quad
    from scipy.special import logsumexp, roots_legendre

    x = _validated_effects(source_effects)
    n = np.asarray(source_n, dtype=float)

    if (n.shape != x.shape or not np.isfinite(n).all()
            or np.any(n < 2) or np.any(n != np.floor(n))):
        raise ValueError(
            "source_n must contain matching integer sample sizes >= 2"
        )
    if not np.isfinite(min_scale) or min_scale <= 0:
        raise ValueError(
            "min_scale must be positive and fixed before fitting"
        )

    root_n = np.sqrt(n)
    t = root_n * x
    log_floor = math.log(min_scale)
    rules = {}

    # Geometric endpoint panels resolve small-scale half-Cauchy tails.
    edges = np.r_[
        0., 10.**np.arange(-6, 0), .5,
        1 - 10.**np.arange(-1, -7, -1), 1.
    ]

    def log_marginals(parameters, order):
        location, log_scale = parameters

        if (not np.isfinite(parameters).all()
                or log_scale < log_floor - 1e-12
                or log_scale > 300):
            return None

        scale = math.exp(log_scale)

        if order not in rules:
            nodes, weights = roots_legendre(order)
            width = np.diff(edges)
            u = (
                (edges[:-1] + edges[1:])[:, None] / 2
                + width[:, None] * nodes / 2
            ).ravel()

            rules[order] = (
                np.tan(math.pi * u / 2),
                np.log((width[:, None] * weights / 2).ravel()),
            )

        lam, log_weights = rules[order]
        a = np.hypot(
            1., root_n[:, None] * scale * lam[None, :]
        )

        # Integrate the normal part of the horseshoe analytically:
        # conditional on lambda, T/a is noncentral t.
        terms = stats.nct.logpdf(
            t[:, None] / a,
            n[:, None] - 1,
            root_n[:, None] * location / a,
        ) - np.log(a)

        if np.isnan(terms).any() or np.isposinf(terms).any():
            return None

        answer = logsumexp(
            terms + log_weights[None, :], axis=1
        )
        return answer if np.isfinite(answer).all() else None

    def objective(parameters, order):
        logs = log_marginals(parameters, order)
        return np.inf if logs is None else -float(logs.sum())

    spread = max(float(np.std(x, ddof=1)), 10 * min_scale)
    starts = [
        np.array([
            location,
            math.log(max(spread * factor, min_scale)),
        ])
        for location in np.unique([np.mean(x), np.median(x)])
        for factor in (.1, .5, 1.)
    ]

    converged = []
    for start in starts:
        result = optimize.minimize(
            objective,
            start,
            args=(12,),
            method="Nelder-Mead",
            bounds=[(None, None), (log_floor, None)],
            options={
                "maxiter": 2000,
                "xatol": 1e-8,
                "fatol": 1e-9,
            },
        )
        if result.success and np.isfinite(result.fun):
            converged.append(result)

    if not converged:
        raise RuntimeError(
            "All horseshoe marginal-likelihood fits failed"
        )

    best = min(converged, key=lambda result: result.fun)

    # Refit at higher resolution and check every source contribution.
    for order in (24, 48, 96):
        refined = optimize.minimize(
            objective,
            best.x,
            args=(order,),
            method="Nelder-Mead",
            bounds=[(None, None), (log_floor, None)],
            options={
                "maxiter": 2000,
                "xatol": 1e-8,
                "fatol": 1e-9,
            },
        )

        if not refined.success or not np.isfinite(refined.fun):
            raise RuntimeError(
                "Horseshoe marginal-likelihood refinement failed"
            )

        best = refined
        logs = log_marginals(best.x, order)
        check_logs = log_marginals(best.x, 2 * order)

        if logs is None or check_logs is None:
            raise FloatingPointError(
                "Nonfinite horseshoe marginal density"
            )

        quadrature_error = float(
            np.max(np.abs(logs - check_logs))
        )
        if quadrature_error < 1e-7:
            break
    else:
        raise FloatingPointError(
            "Horseshoe marginal quadrature did not stabilize"
        )

    location, log_scale = map(float, best.x)
    scale = math.exp(log_scale)

    # Independent adaptive-integration check.
    adaptive_logs = []
    for ti, ni, rn in zip(t, n, root_n):

        def integrand(u):
            if u == 1.:
                return 0.
            a = math.hypot(
                1., rn * scale * math.tan(math.pi * u / 2)
            )
            return float(
                stats.nct.pdf(
                    ti / a, ni - 1, rn * location / a
                ) / a
            )

        result = quad(
            integrand,
            0., 1.,
            points=edges[1:-1],
            epsabs=1e-12,
            epsrel=1e-9,
            limit=250,
            full_output=1,
        )
        value, error = result[:2]

        if (len(result) != 3
                or not np.isfinite([value, error]).all()
                or value <= 0
                or error > 1e-7 * value):
            raise FloatingPointError(
                "Adaptive horseshoe marginal check failed"
            )

        adaptive_logs.append(math.log(value))

    adaptive_error = float(
        np.max(np.abs(check_logs - adaptive_logs))
    )
    if adaptive_error > 1e-6:
        raise FloatingPointError(
            "Horseshoe marginal integration checks disagree"
        )

    at_floor = scale <= min_scale * (1 + 1e-4)
    if at_floor:
        warnings.warn(
            "Horseshoe scale reached its lower bound; report the "
            "boundary fit and check sensitivity to min_scale.",
            RuntimeWarning,
        )

    details = dict(
        actual_method=(
            "constrained marginal MLE with noncentral-t sampling"
        ),
        likelihood_scope=(
            "full real-line latent horseshoe; condition after fitting"
        ),
        input_interpretation=(
            "source t-statistics with their sample sizes"
        ),
        source_weighting=(
            "one marginal likelihood term per source lineage"
        ),
        marginal_nll=-float(np.sum(adaptive_logs)),
        scale_lower_bound=float(min_scale),
        scale_at_lower_bound=at_floor,
        optimizer="multistart Nelder-Mead in location/log-scale",
        optimizer_converged=True,
        starts_attempted=len(starts),
        starts_converged=len(converged),
        global_optimum_certified=False,
        fit_quadrature_order_per_panel=order,
        fit_quadrature_panels=len(edges) - 1,
        marginal_logpdf_quadrature_error=quadrature_error,
        marginal_logpdf_adaptive_error=adaptive_error,
    )

    return Prior("Horseshoe", location, scale), details


def build_priors(source_effects, source_n):
    """Return the seven priors and their source-fitting diagnostics.

    The input is one mean(-Chronos)/sd(Chronos) estimate per source lineage.
    Four families use full-family MLE; Horseshoe uses median/IQR matching.
    Rouder and improper flat are fixed comparisons, never fitted.
    """
    x = _validated_effects(source_effects)
    priors, diagnostics = {}, []
    for family in FAMILIES:
        initializer = (
            fit_prior(x, family) if family != "Horseshoe" else None
        )
        info = dict(
            prior=family,
            requested_method=(
                "four full-family MLE fits; horseshoe marginal MLE"),
            source_lineages=len(x),
            source_weighting='equal lineage weights',
            likelihood_scope='full real-line family; condition on delta>0 after fitting',
            input_interpretation='source effect estimates treated as observations; estimation error not modelled',
        )
        if family != 'Horseshoe':
            prior, numerical = _scipy_mle(x, family, initializer)
            info.update(actual_method='full-family MLE', **numerical)
        else:
            prior, numerical = fit_horseshoe_mml(x, source_n)
            info.update(numerical)

        priors[family] = prior
        info.update(location_positive_delta=prior.location, scale=prior.scale,
                    positive_prior_mass=prior.positive_mass)
        diagnostics.append(info)

    prior = Prior('Cauchy', 0., 1/math.sqrt(2))
    priors[ROUDER_LABEL] = prior
    diagnostics.append(dict(
        prior=ROUDER_LABEL, requested_method='fixed comparison',
        actual_method='fixed proper comparison prior; not fitted',
        source_lineages=0, location_positive_delta=prior.location,
        scale=prior.scale, positive_prior_mass=prior.positive_mass,
        optimizer_converged=None,
    ))
    priors[IMPROPER_LABEL] = FlatPrior()
    diagnostics.append(dict(
        prior=IMPROPER_LABEL, requested_method='fixed comparison',
        actual_method='fixed improper Lebesgue measure on delta>0; not fitted',
        source_lineages=0, location_positive_delta=None, scale=None,
        positive_prior_mass=None, total_mass='infinite',
        optimizer_converged=None,
    ))
    return priors, diagnostics

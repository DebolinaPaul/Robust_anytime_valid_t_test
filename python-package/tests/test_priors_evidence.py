"""Numerical references exported independently from the simulation notebook."""

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import exp1
from scipy.stats import cauchy, laplace, nct, norm, t as student_t

from eville import ev_log_mixture, ev_prior, ev_prior_density


FAMILIES = ["gaussian", "cauchy", "student_t5", "laplace", "horseshoe"]


def test_all_prior_densities_match_notebook(reference):
    for row in reference("all_priors_densities.csv").itertuples():
        prior = ev_prior(row.family, location=row.location, scale=row.scale)
        assert ev_prior_density(row.theta, prior, log=True) == pytest.approx(
            row.log_density, abs=2e-9, rel=0
        )
    assert math.isinf(ev_prior_density(-2, ev_prior("horseshoe", -2, 1.3)))
    assert ev_prior_density(0.7, ev_prior("rouder"), log=True) == pytest.approx(
        cauchy.logpdf(0.7, scale=1 / math.sqrt(2)), abs=1e-14
    )


def test_all_mixture_evidence_matches_notebook(reference):
    for row in reference("all_priors_evidence.csv").itertuples():
        prior = ev_prior(row.family, location=row.location, scale=row.scale)
        assert ev_log_mixture(row.t, row.n, prior) == pytest.approx(
            row.log_evidence, abs=2e-7, rel=0
        )


def test_directional_evidence_matches_existing_cross_language_references(reference):
    data = reference("shifted_evidence.csv")
    # All shifted/directional cases at the smallest and largest n cover the initial
    # calibration distribution as well as evidence after more observations.
    ns = sorted(data.n.unique())
    selected = data[data.n.isin([ns[0], ns[-1]])]
    for row in selected.itertuples():
        kwargs = {"alternative": row.alternative}
        if row.prior != "lebesgue":
            kwargs.update(location=row.location, scale=row.scale)
        prior = ev_prior(row.prior, **kwargs)
        assert ev_log_mixture(row.t, row.n, prior) == pytest.approx(
            row.log_m, abs=3e-7, rel=0
        )


@pytest.mark.parametrize("family", FAMILIES)
def test_one_sided_mixtures_recombine_and_reflect(family):
    full = ev_prior(family, location=-2.3, scale=0.7)
    negative = ev_prior(family, -2.3, 0.7, alternative="less")
    positive = ev_prior(family, -2.3, 0.7, alternative="greater")
    reflected = ev_prior(family, 2.3, 0.7, alternative="greater")
    recombined = np.logaddexp(
        full.log_negative_mass + ev_log_mixture(1.2, 8, negative),
        full.log_positive_mass + ev_log_mixture(1.2, 8, positive),
    )
    assert recombined == pytest.approx(ev_log_mixture(1.2, 8, full), abs=3e-7)
    assert ev_log_mixture(1.2, 8, negative) == pytest.approx(
        ev_log_mixture(-1.2, 8, reflected), abs=3e-7
    )
    assert ev_prior_density(0.4, negative) == 0
    assert ev_prior_density(-0.4, positive) == 0


@pytest.mark.parametrize("family", FAMILIES)
def test_directional_integral_against_noncentral_t_density(family):
    # Independent check: integrate the noncentral-t density ratio directly,
    # rather than evaluating the package's invariant-likelihood quadrature.
    location, scale, n, statistic = -1.2, 0.8, 7, -1.3
    distributions = {"gaussian": norm, "cauchy": cauchy,
                     "student_t5": student_t, "laplace": laplace}

    if family == "horseshoe":
        def density(theta):
            a = 0.5 * ((theta - location) / scale) ** 2
            return math.exp(a) * exp1(a) / (scale * math.pi * math.sqrt(2 * math.pi))
        # Independent latent half-normal representation of the positive tail.
        tail = quad(lambda z: 2 * norm.pdf(z) * math.atan(z / abs(location / scale)) / math.pi,
                    0, np.inf, epsabs=1e-12)[0]
        mass = 1 - tail
    else:
        distribution = distributions[family]
        args = (5,) if family == "student_t5" else ()
        density = lambda theta: distribution.pdf(theta, *args, loc=location, scale=scale)
        mass = distribution.cdf(0, *args, loc=location, scale=scale)

    denominator = student_t.pdf(statistic, n - 1)
    def integrand(theta):
        return nct.pdf(statistic, n - 1, math.sqrt(n) * theta) / denominator * density(theta) / mass

    # The omitted interval theta < -25 has a Gaussian likelihood envelope
    # and is negligible here. Splitting at the horseshoe spike avoids ever
    # evaluating the infinite density at a quadrature node.
    expected = sum(quad(integrand, a, b, epsabs=1e-10, epsrel=1e-9, limit=250)[0]
                   for a, b in [(-25, location), (location, 0)])
    got = ev_log_mixture(statistic, n, ev_prior(family, location, scale, alternative="less"))
    assert got == pytest.approx(math.log(expected), abs=4e-7, rel=0)


def test_half_lebesgue_recombination_and_n2_formula():
    full = ev_prior("lebesgue")
    pos = ev_prior("lebesgue", alternative="greater")
    neg = ev_prior("lebesgue", alternative="less")
    for n in (2, 7, 50):
        for statistic in (-50.0, -1.0, 0.0, 2.0, 50.0):
            lp = ev_log_mixture(statistic, n, pos)
            ln = ev_log_mixture(statistic, n, neg)
            assert np.logaddexp(lp, ln) == pytest.approx(
                ev_log_mixture(statistic, n, full), abs=2e-8
            )
            if n == 2:
                expected = 0.5 * math.log(math.pi) - math.log(4) + np.logaddexp(0, 2 * np.arcsinh(statistic))
                assert lp == pytest.approx(expected, abs=1e-12)


def test_horseshoe_extreme_shift_tail_normalization():
    # Regression for the historical shifted-horseshoe survival-normalization
    # failure: P(lambda*Z > x) ~ sqrt(2)/(pi**1.5*x).
    prior = ev_prior("horseshoe", location=-1e8, scale=1.0, alternative="greater")
    expected_log_mass = math.log(math.sqrt(2) / math.pi ** 1.5) - math.log(1e8)
    assert prior.log_positive_mass == pytest.approx(expected_log_mass, abs=1e-8)
    assert math.isfinite(ev_prior_density(1.0, prior, log=True))


def test_horseshoe_shifted_directional_spike_regression():
    # Independently integrated reference for a spike inside the integration
    # range; ordinary quadrature can converge slowly at the log singularity.
    prior = ev_prior("horseshoe", location=-3, alternative="less")
    assert ev_log_mixture(-1.8708684117893888, 2, prior) == pytest.approx(
        0.18407380204843393, abs=2e-9, rel=0
    )


@pytest.mark.parametrize("kwargs", [
    {"family": "unknown"}, {"family": "gaussian", "scale": 0},
    {"family": "cauchy", "location": np.nan},
    {"family": "rouder", "location": 1},
    {"family": "lebesgue", "scale": 1},
    {"family": "gaussian", "alternative": "left"},
])
def test_invalid_prior_inputs_rejected(kwargs):
    with pytest.raises((ValueError, TypeError)):
        ev_prior(**kwargs)

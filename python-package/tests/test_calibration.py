"""Calibration checks distinguish numerical brackets from sampling error."""

import math

import numpy as np
import pytest
from scipy.integrate import quad

from eville import ev_calibration, ev_log_mixture, ev_null_sample, ev_prior


def test_bundled_null_sample_matches_notebook_and_does_not_touch_rng():
    np.random.seed(17)
    state = np.random.get_state()
    sample = ev_null_sample()
    after = np.random.get_state()
    assert len(sample) == 100000
    np.testing.assert_allclose(sample[:5], [49.61297690645055, 9.423664797088913,
        -5.850566194162971, -0.36620651237260304, -1.537296283600286], atol=1e-14, rtol=0)
    np.testing.assert_array_equal(sample[:257], ev_null_sample(257))
    assert state[0] == after[0]
    np.testing.assert_array_equal(state[1], after[1])
    assert state[2:] == after[2:]
    # The caller cannot corrupt the cache by editing the returned array.
    try:
        sample[0] = 0
    except ValueError:
        pass
    assert ev_null_sample(2)[0] == pytest.approx(49.61297690645055)


def test_empirical_calibration_bounds_contain_exact_sample_average(reference):
    ref = reference("all_priors_mc_direct.csv")
    for family, rows in ref.groupby("family", sort=False):
        first = rows.iloc[0]
        prior = ev_prior(family, first.location, first.scale)
        calibration = ev_calibration(prior, null_t=ev_null_sample(257), tolerance=2e-7)
        for row in rows.itertuples():
            result = calibration.evaluate(row.log_c, relative_tolerance=None)
            assert result["B"] == 257
            assert result["mc_lower"] <= row.empirical_mean + 5e-9
            assert result["mc_upper"] >= row.empirical_mean - 5e-9
            assert result["mc_gap"] <= 2.01e-7
            assert result["mc_estimate"] == pytest.approx((result["mc_lower"] + result["mc_upper"]) / 2, abs=1e-12)
            assert result["p_eville"] == pytest.approx(min(result["mc_upper"], row.ville_bound), abs=1e-12)
            assert result["mcse"] == pytest.approx(row.direct_mcse, abs=2e-6)


def test_shared_100000_draw_calibration_matches_notebook(reference):
    ref = reference("all_priors_mc_100000.csv")
    for family, rows in ref.groupby("family", sort=False):
        first = rows.iloc[0]
        calibration = ev_calibration(ev_prior(family, first.location, first.scale), tolerance=first.tolerance)
        for row in rows.itertuples():
            result = calibration.evaluate(row.log_c, relative_tolerance=None)
            assert result["B"] == 100000
            assert result["mc_lower"] <= row.upper + 2e-8
            assert result["mc_upper"] >= row.lower - 2e-8
            assert result["p_eville"] == pytest.approx(row.used, abs=2.1e-5)
            assert result["mcse"] == pytest.approx(row.mcse, abs=2e-6)
            assert result["p_eville"] <= row.ville_bound + 1e-14


@pytest.mark.parametrize("alternative", ["two.sided", "less", "greater"])
def test_lebesgue_analytic_calibration_matches_null_angle_integral(alternative):
    prior = ev_prior("lebesgue", alternative=alternative)
    calibration = ev_calibration(prior)
    assert calibration.method == "analytic"
    for c in (0.3, 2.0, 20.0, 200.0):
        # Under H0, T2 is Cauchy; its angle is uniform on (-pi/2, pi/2).
        def integrand(angle):
            logm = ev_log_mixture(math.tan(angle), 2, prior)
            return math.exp(min(0, logm - math.log(c))) / math.pi
        expected = quad(integrand, -math.pi / 2, math.pi / 2,
                        epsabs=1e-10, points=[0], limit=300)[0]
        assert calibration.evaluate(math.log(c))["p_eville"] == pytest.approx(expected, abs=2e-8)


def test_centered_gaussian_population_calibration_is_analytic():
    scale = 1.3
    calibration = ev_calibration(ev_prior("gaussian", scale=scale), method="deterministic")
    b = 2 * scale ** 2
    a = math.sqrt(1 + b)
    for c in (0.2, 0.8, 1.1, 1.7, 20):
        if c <= 1 / a:
            expected = 1.0
        elif c >= a:
            expected = 1 / c
        else:
            cosine = math.sqrt((a / c - 1) / b)
            sine = math.sqrt(1 - cosine ** 2)
            expected = 2 / math.pi * (math.atan2(sine, a * cosine) / c + math.asin(cosine))
        assert calibration.evaluate(math.log(c))["p_eville"] == pytest.approx(expected, abs=2e-8)


@pytest.mark.parametrize("family", ["gaussian", "cauchy"])
@pytest.mark.parametrize("alternative", ["less", "greater"])
def test_general_population_calibration_brackets_independent_references(reference, family, alternative):
    ref = reference("shifted_calibration.csv")
    rows = ref[(ref.prior == family) & (ref.location == -3) & (ref.alternative == alternative)]
    calibration = ev_calibration(
        ev_prior(family, location=-3, scale=1.4, alternative=alternative),
        method="deterministic", tolerance=1e-5,
    )
    for row in rows.itertuples():
        result = calibration.evaluate(row.log_c)
        assert result["lower"] <= row.g + 3e-7
        assert result["p_eville"] >= row.g - 3e-7
        assert result["gap"] <= 1.01e-5


@pytest.mark.parametrize("alternative", ["less", "greater"])
def test_one_sided_mc_brackets_direct_empirical_mean(alternative):
    prior = ev_prior("gaussian", location=-1.5, scale=0.9, alternative=alternative)
    null = ev_null_sample(31)
    logs = np.array([ev_log_mixture(t, 2, prior) for t in null])
    calibration = ev_calibration(prior, null_t=null, tolerance=1e-7)
    for c in (0.3, 2.0, 20.0):
        direct = np.exp(np.minimum(0, logs - math.log(c))).mean()
        result = calibration.evaluate(math.log(c), relative_tolerance=None)
        assert result["mc_lower"] <= direct + 5e-9
        assert result["mc_upper"] >= direct - 5e-9


def test_calibration_invalid_inputs_rejected():
    with pytest.raises((ValueError, TypeError)):
        ev_calibration(ev_prior(), null_t=[0, np.nan])
    with pytest.raises((ValueError, TypeError)):
        ev_calibration(ev_prior(), method="deterministic", null_t=[-1, 1])
    with pytest.raises((ValueError, TypeError)):
        ev_calibration(ev_prior("lebesgue"), null_t=[-1, 1])
    with pytest.raises((ValueError, TypeError)):
        ev_null_sample(1)

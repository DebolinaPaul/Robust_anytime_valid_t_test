"""Public one-sample test, result summaries, plotting and transformations."""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from eville import ev_calibration, ev_log_evidence, ev_null_sample, ev_prior, ev_ttest


OBSERVATIONS = np.array([-0.4, 0.9, 0.1, 1.3, 0.4, 1.2, 0.6, 1.4])


@pytest.mark.parametrize("family", ["gaussian", "cauchy", "student_t5", "laplace", "horseshoe"])
def test_paths_preserve_scale_invariance_and_same_prior_ordering(family):
    prior = ev_prior(family, location=2)
    calibration = ev_calibration(prior, null_t=ev_null_sample(257))
    fit = ev_ttest(OBSERVATIONS, prior=prior, calibration=calibration)
    scaled = ev_ttest(1e60 * OBSERVATIONS, prior=prior, calibration=calibration)
    path = fit.to_dataframe()
    np.testing.assert_array_equal(path.n, np.arange(2, len(OBSERVATIONS) + 1))
    np.testing.assert_allclose(path.log_evidence, scaled.path.log_evidence, atol=2e-8, rtol=0)
    assert (path.p_eville <= path.p_ville + 1e-14).all()
    assert (np.diff(path.p_eville) <= 1e-14).all()
    assert (np.diff(path.p_ville) <= 1e-14).all()
    for method, column in [("eVille", "p_eville"), ("Ville", "p_ville")]:
        crossed = path.loc[path[column] <= fit.alpha, "n"]
        assert fit.stopping_times[method] == (None if crossed.empty else int(crossed.iloc[0]))
    assert ev_log_evidence(OBSERVATIONS, prior) == pytest.approx(path.log_evidence.iloc[-1], abs=1e-10)


def test_null_mean_shift_and_directional_reflection():
    prior = ev_prior("gaussian", location=-1, alternative="less")
    reflected = ev_prior("gaussian", location=1, alternative="greater")
    null = ev_null_sample(31)
    cal = ev_calibration(prior, null_t=null)
    # Reflect the calibration sample as well; empirical MC is not exactly
    # reflection invariant on a nonsymmetric finite sample.
    cal_reflected = ev_calibration(reflected, null_t=-null)
    first = ev_ttest(OBSERVATIONS, prior=prior, calibration=cal)
    shifted = ev_ttest(OBSERVATIONS + 3, prior=prior, mu=3, calibration=cal)
    mirror = ev_ttest(-OBSERVATIONS, prior=reflected, calibration=cal_reflected)
    np.testing.assert_allclose(first.path.log_evidence, shifted.path.log_evidence, atol=2e-8)
    np.testing.assert_allclose(first.path.log_evidence, mirror.path.log_evidence, atol=2e-8)
    np.testing.assert_allclose(first.path.p_eville, mirror.path.p_eville, atol=3e-5)


def test_summaries_and_plots_handle_unavailable_ville(tmp_path):
    for family in ("gaussian", "lebesgue"):
        prior = ev_prior(family)
        fit = ev_ttest(OBSERVATIONS, prior=prior, calibration_method="deterministic")
        summary = fit.summary()
        assert set(summary["method"]) == {"Ville", "eVille"}
        if family == "lebesgue":
            assert fit.path.p_ville.isna().all()
            assert fit.stopping_times["Ville"] is None
            assert "unavailable" in str(summary.loc[summary.method == "Ville", "decision"].iloc[0]).lower()
        axes = fit.plot()
        axes.figure.savefig(tmp_path / f"{family}.pdf")
        assert (tmp_path / f"{family}.pdf").stat().st_size > 1000
        plt.close(axes.figure)


@pytest.mark.parametrize("x", [[1], [1, 1, 2], [1, np.nan], [1, np.inf], [[1, 2], [3, 4]]])
def test_invalid_observations_rejected(x):
    with pytest.raises((ValueError, TypeError)):
        ev_ttest(x)


def test_reject_mismatched_calibration_and_invalid_level():
    cal = ev_calibration(ev_prior("gaussian"), method="deterministic")
    with pytest.raises((ValueError, TypeError)):
        ev_ttest(OBSERVATIONS, prior=ev_prior("gaussian", location=1), calibration=cal)
    with pytest.raises((ValueError, TypeError)):
        ev_ttest(OBSERVATIONS, n0=5)
    with pytest.raises((ValueError, TypeError)):
        ev_ttest(OBSERVATIONS, alpha=0)
    with pytest.raises((ValueError, TypeError)):
        ev_ttest(OBSERVATIONS, prior=ev_prior("gaussian"), calibration=cal, null_t=[-1, 1])

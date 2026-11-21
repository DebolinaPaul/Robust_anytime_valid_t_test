"""Compare five fixed priors on one reproducible synthetic sequence."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from eville import ev_calibration, ev_null_sample, ev_prior, ev_ttest


def main():
    output = Path("example_outputs")
    output.mkdir(exist_ok=True)

    # This target seed and the default prior scales match the notebook.
    # All priors below are specified in advance of observing this sequence.
    x = np.random.default_rng(20260509).normal(0.5, 1, 50)
    null_t = ev_null_sample()
    families = ["gaussian", "cauchy", "student_t5", "laplace", "horseshoe"]
    summaries = []

    for family in families:
        prior = ev_prior(family, location=3)
        calibration = ev_calibration(prior, null_t=null_t)
        fit = ev_ttest(x, prior=prior, calibration=calibration, alpha=0.05)
        summary = fit.summary().copy()
        summary.insert(0, "prior_family", family)
        summaries.append(summary)
        fit.to_dataframe().to_csv(output / f"{family}_trajectory.csv", index=False)
        ax = fit.plot()
        ax.figure.savefig(output / f"{family}_trajectory.pdf", bbox_inches="tight")
        plt.close(ax.figure)

    comparison = pd.concat(summaries, ignore_index=True)
    comparison.to_csv(output / "prior_comparison.csv", index=False)
    pd.DataFrame({"n": np.arange(1, len(x) + 1), "x": x}).to_csv(
        output / "target_observations.csv", index=False
    )
    print(comparison.to_string(index=False))
    print("eVille results use finite Monte Carlo calibration.")
    print(f"Results saved in {output.resolve()}")


if __name__ == "__main__":
    main()

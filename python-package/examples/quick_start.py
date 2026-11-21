"""Run a one-sided one-sample t-test and export its results."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from eville import ev_prior, ev_ttest


def main():
    output = Path("example_outputs")
    output.mkdir(exist_ok=True)

    # Synthetic observations in their sampling order. Specify the prior
    # independently of these target observations.
    x = [-0.3, -1.1, -0.7, -0.5, -1.4, -0.8, -0.2, -1.0]
    prior = ev_prior("gaussian", location=-4, scale=1, alternative="less")
    fit = ev_ttest(x, prior=prior, alpha=0.05)

    print(fit.summary().to_string(index=False))
    fit.summary().to_csv(output / "quick_start_summary.csv", index=False)
    fit.to_dataframe().to_csv(output / "quick_start_trajectory.csv", index=False)
    ax = fit.plot()
    ax.figure.savefig(output / "quick_start.pdf", bbox_inches="tight")
    plt.close(ax.figure)
    print(f"Results saved in {output.resolve()}")


if __name__ == "__main__":
    main()

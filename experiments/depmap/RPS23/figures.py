#!/usr/bin/env python3
"""Regenerate the four paper figures from consolidated result CSVs.

Run ``python figures.py --results results`` after run_analysis.py.
This module reads saved results; it does not refit priors or rerun tests.
Random-order summaries describe order sensitivity of the observed samples.
"""
from __future__ import annotations
import pandas as pd
import numpy as np
from matplotlib.ticker import MaxNLocator, PercentFormatter
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")


PROPER_PRIORS = ["Gaussian", "Cauchy", "Laplace", "Student-t5",
                 "Horseshoe", "Rouder half-Cauchy"]
FLAT = "Improper flat"
PRIORS = PROPER_PRIORS + [FLAT]
PRIOR_LABELS = {"Gaussian": "Gaussian", "Cauchy": "Cauchy",
                "Laplace": "Laplace", "Student-t5": r"Student-$t_5$",
                "Horseshoe": "Horseshoe", "Rouder half-Cauchy": "Rouder",
                FLAT: "Improper flat"}
EVILLE = "#16775B"
VILLE = "#7B3E8D"
TARGET = "#C66A26"
INK = "#24333C"


def configure_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8,
        "axes.titlesize": 9, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 8, "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#7B858C", "axes.linewidth": .6,
        "text.color": INK, "axes.labelcolor": INK,
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })


def _save(fig, out, name):
    path = out / name
    temporary = out / (path.stem + ".tmp.pdf")
    fig.savefig(temporary, format="pdf")
    plt.close(fig)
    temporary.replace(path)
    return path


def _boolean(values):
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).to_numpy(bool)
    mapped = values.astype(str).str.lower().map({"true": True, "false": False,
                                                "1": True, "0": False})
    if mapped.isna().any():
        raise ValueError("An eligibility column contains an invalid boolean")
    return mapped.to_numpy(bool)


def _lineages(summary):
    records = summary[["target_lineage",
                       "target_n", "cohort"]].drop_duplicates()
    if records.target_lineage.duplicated().any():
        raise ValueError("Lineage sizes or cohort labels are inconsistent")
    records = records.assign(
        _group=records.cohort.map({"main": 0, "small": 1})
    )
    if records._group.isna().any():
        raise ValueError(
            "Expected main/small cohort labels for tested lineages")
    return records.sort_values(["_group", "target_lineage"]).reset_index(drop=True)


def paired_savings_all_lineages(summary, settings, out):
    """Show all paired outcomes; nonrejections never become zero-saving ties."""
    rows = _lineages(summary)
    names = rows.target_lineage.tolist()
    labels = [
        f"{r.target_lineage} ({int(r.target_n)})" for r in rows.itertuples()]
    nmain = int(rows.cohort.eq("main").sum())
    categories = [
        ("outcome_same_rejection_fraction", "Same rejection time", "#D9DEE2", None),
        ("outcome_save1_fraction", "Save 1 observation", "#A9D6C4", None),
        ("outcome_save2_fraction", "Save 2 observations", "#54A687", None),
        ("outcome_save3plus_fraction", "Save 3+ observations", "#17684D", None),
        ("eville_only_rejects_fraction", "Only eVille rejects", "#E7A15B", None),
        ("neither_rejects_fraction", "Neither rejects", "white", "////"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(
        7.5, 9.1), sharex=True, sharey=True)
    fig.subplots_adjust(left=.245, right=.972, top=.925, bottom=.145,
                        wspace=.17, hspace=.22)
    y = np.arange(len(names))
    for index, (prior, ax) in enumerate(zip(PROPER_PRIORS, axes.flat)):
        selected = summary.loc[summary.prior.eq(
            prior)].set_index("target_lineage")
        if selected.index.duplicated().any():
            raise ValueError(f"Duplicate lineage rows for {prior}")
        matrix = selected.reindex(names)[[c[0]
                                          for c in categories]].to_numpy(float)
        if (not np.isfinite(matrix).all() or np.any(matrix < -1e-12)
                or not np.allclose(matrix.sum(axis=1), 1., atol=1e-10)):
            raise ValueError(
                f"Invalid paired outcome probabilities for {prior}")
        left = np.zeros(len(names))
        for j, (_, _, color, hatch) in enumerate(categories):
            ax.barh(y, matrix[:, j], left=left, height=.72, color=color,
                    edgecolor="#7B858C" if hatch else "white",
                    linewidth=.25, hatch=hatch)
            left += matrix[:, j]
        if 0 < nmain < len(names):
            ax.axhline(nmain-.5, color="#66747D", ls=(0, (3, 3)), lw=.7)
        ax.set_xlim(0, 1)
        ax.set_ylim(len(names)-.55, -.65)
        ax.set_xticks([0, .5, 1])
        ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.set_yticks(y)
        ax.tick_params(axis="y", length=0, labelleft=index % 3 == 0)
        if index % 3 == 0:
            ax.set_yticklabels(labels, fontsize=7.2)
        ax.tick_params(axis="x", length=2.5)
        ax.spines["left"].set_visible(False)
        ax.set_title(
            f"({chr(97+index)}) {PRIOR_LABELS[prior]}", loc="left", pad=7)
        if index >= 3:
            ax.set_xlabel("Paired observation orders")
    fig.suptitle("Paired stopping outcomes across held-out lineages",
                 x=.025, y=.98, ha="left", fontsize=11)
    fig.text(.025, .95,
             f"{settings['gene']}; recorded Sex = {settings['sex']}; "
             f"{settings['permutations']} random orders per lineage", fontsize=8)
    handles = [Patch(facecolor=color, edgecolor="#7B858C" if hatch else "white",
                     linewidth=.4, hatch=hatch, label=label)
               for _, label, color, hatch in categories]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .056),
               ncol=3, frameon=False, fontsize=7.5, columnspacing=1.3,
               handlelength=1.7)
    fig.text(.025, .016,
             "Numbers in parentheses are target sample sizes. Below dashed line: small-target sensitivity analysis.\n"
             "Saving categories require both methods to reject. Horseshoe uses marginal maximum likelihood; improper flat has no Ville comparison.",
             fontsize=7, linespacing=1.45)
    return _save(fig, out, "paired_savings_all_lineages.pdf")


def conflict_and_flat_baseline(
    summary, settings, out, fits, *,
    lineage=None,
    filename="conflict_and_flat_baseline.pdf",
):
    standardize = True  # False: use |theta_hat - prior location|.
    fitted = [p for p in PROPER_PRIORS if p != "Rouder half-Cauchy"]
    main = summary.loc[
        summary.cohort.eq("main") & summary.target_n.ge(10)
    ].copy()

    data = main.loc[main.prior.isin(fitted)].merge(
        fits[["target_lineage", "prior",
              "full_prior_location_raw_theta", "scale"]],
        on=["target_lineage", "prior"],
        validate="one_to_one",
    )
    data["x"] = (
        data.theta_hat - data.full_prior_location_raw_theta
    ).abs()
    if standardize:
        data["x"] /= data.scale

    rouder = main.loc[
        main.prior.eq("Rouder half-Cauchy")
    ].set_index("target_lineage")
    flat = main.loc[
        main.prior.eq(FLAT)
    ].set_index("target_lineage")

    same_rouder = np.allclose(
        rouder.eville_mean_capped_at_N,
        rouder.ville_mean_capped_at_N,
        rtol=0,
        atol=1e-12,
    )
    values = np.r_[
        main.eville_mean_capped_at_N,
        main.loc[main.prior.ne(FLAT), "ville_mean_capped_at_N"],
    ]
    if not len(data) or not np.isfinite(values).all():
        raise ValueError("Missing or nonfinite main-cohort results")
    span = max(float(np.ptp(values)), .5)
    # Filter after calculating the common vertical-axis limits.
    if lineage is not None:
        data = data.loc[data.target_lineage.eq(lineage)].copy()
        if data.empty:
            raise ValueError(
                f"No eligible results for lineage: {lineage}"
            )

    fig, axes = plt.subplots(2, 3, figsize=(7.5, 5.7), sharey=True)
    fig.subplots_adjust(
        left=.10, right=.98, top=.88, bottom=.20,
        wspace=.18, hspace=.38,
    )

    for i, (prior, ax) in enumerate(zip(fitted, axes.flat)):
        d = data.loc[data.prior.eq(prior)].sort_values(
            ["x", "target_lineage"], kind="stable"
        )
        x = d.x.to_numpy(float)
        names = d.target_lineage.tolist()
        r, f = rouder.loc[names], flat.loc[names]

        series = [
            (d.eville_mean_capped_at_N, EVILLE, "-", "o",
             "Fitted prior: eVille"),
            (d.ville_mean_capped_at_N, VILLE, "--", "x",
             "Fitted prior: Ville"),
            (r.eville_mean_capped_at_N, "#C66A26", "-.", "s",
             "Rouder: eVille" if same_rouder
             else "Rouder: eVille"),
            (f.eville_mean_capped_at_N, "#66747D", ":", "^",
             "Lebesgue: eVille"),
        ]
        if not same_rouder:
            series.append(
                (r.ville_mean_capped_at_N, "#9B5B32", "--", "D",
                 "Rouder: Ville")
            )

        for y, color, style, marker, label in series:
            y = np.asarray(y, dtype=float)
            if not np.isfinite(x).all() or not np.isfinite(y).all():
                raise ValueError(f"Missing plotting values for {prior}")
            ax.plot(
                x, y, color=color, ls=style, marker=marker,
                lw=.9, ms=3.3, alpha=.85, label=label,
            )

        ax.set_ylim(values.min()-.08*span, values.max()+.16*span)
        ax.margins(x=.07)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.grid(axis="y", color="#E9ECEF", lw=.5)
        ax.set_title(
            f"({chr(97+i)}) {PRIOR_LABELS[prior]}", loc="left"
        )

        example = d.loc[
            d.target_lineage.eq(settings["example_target"])
        ]
        if len(example) == 1:
            row = example.iloc[0]
            ax.annotate(
                settings["example_target"],
                (row.x, row.eville_mean_capped_at_N),
                xytext=(-5, 8), textcoords="offset points",
                ha="right", fontsize=6.8,
            )

    key = axes.flat[-1]
    key.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    key.legend(
        handles, labels, loc="center", frameon=False,
        fontsize=7.5, handlelength=3,
    )

    title = "Stopping time and prior discrepancy"
    if lineage is not None:
        title += f": {lineage}"

    fig.suptitle(
        title, x=.025, y=.98, ha="left", fontsize=11
    )

    subtitle = (
        "Main cohort" if lineage is None
        else f"Target: {lineage}"
    )
    fig.text(
        .025, .935,
        f"{subtitle}; {settings['permutations']} paired random orders",
        fontsize=8,
    )
    xlabel = (
        r"$|\widehat\theta_{\mathrm{target}}-m_{\mathrm{prior}}|"
        r"/s_{\mathrm{prior}}$"
        if standardize else
        r"$|\widehat\theta_{\mathrm{target}}-m_{\mathrm{prior}}|$"
    )
    fig.text(
        .54, .12, "Prior discrepancy: " + xlabel,
        ha="center", fontsize=8,
    )
    note = (
        "Each point is a lineage; joining lines do not represent "
        "a controlled conflict experiment."
        if lineage is None else
        "One observed discrepancy per prior; permutations vary order, "
        "not the full-sample effect estimate."
    )

    fig.text(
        .025, .025,
        note
        + "\nRouder and Lebesgue use the same target as the fitted prior. "
        "Horseshoe uses marginal maximum likelihood.",
        fontsize=6.8,
        linespacing=1.5,
    )

    return _save(fig, out, filename)


def example_forest_and_evidence(cohort, paths, fits, settings, out):
    target = settings["example_target"]
    target_rows = cohort.loc[cohort.target_lineage.eq(target)]
    if len(target_rows) != 1:
        raise ValueError(
            "Example target must appear exactly once in cohort.csv")
    source = cohort.loc[_boolean(
        cohort.source_eligible) & cohort.target_lineage.ne(target)]
    source = source.sort_values("theta_hat")
    forest = pd.concat([source, target_rows], ignore_index=True)
    theta = forest.theta_hat.to_numpy(float)
    lower = forest.theta_ci95_low.to_numpy(float)
    upper = forest.theta_ci95_high.to_numpy(float)
    if not np.isfinite(np.r_[theta, lower, upper]).all():
        raise ValueError(
            "Forest plot requires finite effect estimates and intervals")
    gaussian = paths.loc[paths.prior.eq("Gaussian")].sort_values("n")
    thresholds = fits.loc[fits.target_lineage.eq(
        target) & fits.prior.eq("Gaussian")]
    if gaussian.empty or len(thresholds) != 1:
        raise ValueError("Missing or duplicated Gaussian example results")
    log_e = float(thresholds.iloc[0].log_eville_threshold)
    log_v = math.log(float(thresholds.iloc[0].ville_threshold))
    ns = gaussian.n.to_numpy(int)
    logs = gaussian.log_running_evidence.to_numpy(float)
    e_hits, v_hits = ns[logs >= log_e], ns[logs >= log_v]
    e_stop = int(e_hits[0]) if len(e_hits) else None
    v_stop = int(v_hits[0]) if len(v_hits) else None

    fig = plt.figure(figsize=(7.5, 5.35))
    grid = fig.add_gridspec(2, 2, left=.245, right=.985, top=.86, bottom=.20,
                            width_ratios=(1.12, 1.), wspace=.33, hspace=.50)
    ax_f = fig.add_subplot(grid[:, 0])
    ax_full, ax_zoom = fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, 1])
    yy = np.arange(len(forest))
    # Draw endpoints directly: intervals need not be symmetric about the estimate.
    ax_f.hlines(yy[:-1], lower[:-1], upper[:-1], color="#718999", lw=.85)
    ax_f.scatter(theta[:-1], yy[:-1], color="#496E84", s=13, zorder=3)
    ax_f.hlines(yy[-1], lower[-1], upper[-1], color=TARGET, lw=1.3)
    ax_f.scatter(theta[-1], yy[-1], color=TARGET, marker="D", s=27, zorder=4)
    ax_f.axhline(yy[-1]-.5, color="#CFD5D9", lw=.7)
    ax_f.axvline(0., color="#9EA8AE", ls=":", lw=.8)
    labels = [
        f"{r.target_lineage} ({int(r.target_n)})" for r in forest.itertuples()]
    labels[-1] += " *"
    ax_f.set_yticks(yy, labels, fontsize=7.2)
    ax_f.set_ylim(len(forest)-.4, -.6)
    ax_f.tick_params(axis="y", length=0)
    ax_f.spines["left"].set_visible(False)
    ax_f.set_xlim(min(float(lower.min()), 0.)-.4,
                  max(0., float(upper.max()))+.3)
    ax_f.xaxis.set_major_locator(MaxNLocator(4))
    ax_f.set_xlabel(r"Standardized mean Chronos score, $\theta$")
    ax_f.set_title("(a) Source and target effects", loc="left", pad=9)

    for ax in (ax_full, ax_zoom):
        ax.step(ns, logs, where="post", color=INK, lw=1.25, zorder=3)
        ax.axhline(log_e, color=EVILLE, lw=1.1)
        ax.axhline(log_v, color=VILLE, lw=1.1, ls="--")
        ax.set_ylabel(r"$\log\max_{2\leq k\leq n}M_k$")
        ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.grid(axis="y", color="#E9ECEF", lw=.5)
    ax_full.set_xlim(ns.min(), ns.max())
    ax_full.set_title("(b) Gaussian evidence", loc="left", pad=9)
    ax_full.set_xlabel("Observations, n")
    margin = max(.12, .55*abs(log_v-log_e))
    ax_zoom.set_ylim(min(log_e, log_v)-margin, max(log_e, log_v)+margin)
    last_crossing = max([v for v in (e_stop, v_stop)
                        if v is not None], default=int(ns.max()))
    ax_zoom.set_xlim(ns.min(), min(int(ns.max()), max(4, last_crossing+1)))
    for stop, color, label in [(e_stop, EVILLE, "eVille"), (v_stop, VILLE, "Ville")]:
        if stop is not None:
            ax_zoom.axvline(stop, color=color, lw=.75, ls=":")
            ax_zoom.text(stop+.06, ax_zoom.get_ylim()[0]+.04*(ax_zoom.get_ylim()[1]-ax_zoom.get_ylim()[0]),
                         f"{label}: {stop}", fontsize=6.8, rotation=90,
                         color=color, va="bottom", ha="left")
    ax_zoom.set_title("(c) First-crossing region", loc="left", pad=9)
    ax_zoom.set_xlabel("Observations, n")
    fig.suptitle(f"{settings['gene']}: sources and held-out {target}",
                 x=.025, y=.98, ha="left", fontsize=11)
    fig.text(.025, .93, "Ascending ModelID order; the two calibrations use the same evidence process",
             fontsize=8)
    handles = [Line2D([], [], color=INK, lw=1.25, label="Running log evidence"),
               Line2D([], [], color=EVILLE, lw=1.1, label="eVille threshold"),
               Line2D([], [], color=VILLE, lw=1.1, ls="--", label="Ville threshold")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .065),
               ncol=3, frameon=False, fontsize=7.5)
    fig.text(.025, .018,
             "*Held-out target. Forest intervals are marginal fixed-sample 95% Gaussian-model intervals.\n"
             "This fixed order is illustrative; paired random-order results are reported separately.",
             fontsize=7, linespacing=1.4)
    return _save(fig, out, "example_forest_and_evidence.pdf")


def example_running_pvalues(paths, settings, out):
    alpha = float(settings["alpha"])
    positive = paths[["p_eville", "p_ville"]].to_numpy(float).ravel()
    positive = positive[np.isfinite(positive) & (positive > 0)]
    if not len(positive):
        raise ValueError("No finite positive p-values in the example paths")
    floor = max(1e-300, 10. ** math.floor(math.log10(float(positive.min()))))
    fig, axes = plt.subplots(2, 4, figsize=(
        7.5, 4.75), sharex=True, sharey=True)
    fig.subplots_adjust(left=.085, right=.985, top=.83, bottom=.15,
                        wspace=.12, hspace=.35)
    N = int(paths.n.max())
    for index, (prior, ax) in enumerate(zip(PRIORS, axes.flat)):
        data = paths.loc[paths.prior.eq(prior)].sort_values("n")
        if data.empty:
            raise ValueError(f"Missing example path for {prior}")
        ax.step(data.n, np.maximum(data.p_eville, floor), where="post",
                color=EVILLE, lw=1.3)
        if prior != FLAT:
            ax.step(data.n, np.maximum(data.p_ville, floor), where="post",
                    color=VILLE, lw=1.2, ls="--")
        elif data.p_ville.notna().any():
            raise ValueError("Improper flat must not contain Ville p-values")
        ax.axhline(alpha, color="#849099", lw=.7, ls=":")
        ax.set_yscale("log")
        ax.set_ylim(floor, 1.15)
        ax.set_xlim(int(data.n.min()), N)
        ax.set_xticks(sorted(set([2, max(2, N//2), N])))
        ax.set_title(
            f"({chr(97+index)}) {PRIOR_LABELS[prior]}", loc="left", pad=7)
        ax.grid(axis="y", which="major", color="#E9ECEF", lw=.5)
        if index % 4 == 0:
            ax.set_ylabel("Anytime-valid p-value")
        if index >= 4:
            ax.set_xlabel("Observations, n")
    key = axes.flat[-1]
    key.axis("off")
    handles = [Line2D([], [], color=EVILLE, lw=1.3, label="eVille"),
               Line2D([], [], color=VILLE, lw=1.2, ls="--", label="Ville"),
               Line2D([], [], color="#849099", lw=.8, ls=":", label=rf"$\alpha={alpha:g}$")]
    key.legend(handles=handles, loc="upper left", frameon=False, fontsize=8)
    key.text(.04, .35, "Improper flat:\neVille only", transform=key.transAxes,
             fontsize=7.5, va="top", linespacing=1.5)
    fig.suptitle(f"{settings['gene']}, held-out {settings['example_target']}: running p-values",
                 x=.025, y=.98, ha="left", fontsize=11)
    fig.text(.025, .925,
             f"Recorded Sex = {settings['sex']}; ascending ModelID; monitoring from n = {settings['n0']}",
             fontsize=8)
    fig.text(.025, .025,
             "Horseshoe uses marginal maximum likelihood; Gaussian, Cauchy, Laplace and Student-t use source-only full-family MLE.\n"
             "Rouder and improper flat are fixed comparison measures. The same target order is used in every panel.",
             fontsize=7, linespacing=1.45)
    return _save(fig, out, "example_running_pvalues.pdf")


def raw_stopping_times_all_lineages(summary, settings, out):
    """Unadjusted mean first rejection times across random orders."""
    s = summary.loc[summary.target_n.ge(10)].copy()
    if s.empty or s.duplicated(["target_lineage", "prior"]).any():
        raise ValueError(
            "Need nonempty results with one row per lineage/prior"
        )

    fitted = [
        p for p in PROPER_PRIORS if p != "Rouder half-Cauchy"
    ]
    rows = s[["target_lineage", "target_n"]].drop_duplicates()
    rows = rows.sort_values("target_lineage")
    names = rows.target_lineage.tolist()
    labels = [
        f"{r.target_lineage} ({int(r.target_n)})"
        for r in rows.itertuples()
    ]
    y = np.arange(len(names))

    e_col = "eville_mean_first_rejection_time_if_all_orders_reject"
    v_col = "ville_mean_first_rejection_time_if_all_orders_reject"

    def get_values(prior, column):
        d = s.loc[s.prior.eq(prior)].set_index("target_lineage")
        values = d.reindex(names)[column].to_numpy(float)
        if not np.isfinite(values).all():
            raise ValueError(
                f"Missing mean rejection times for {prior}. "
                "This figure requires every random order to reject; "
                "nonrejections cannot be treated as observed "
                "rejection times."
            )
        return values

    fitted_values = {
        p: (get_values(p, e_col), get_values(p, v_col))
        for p in fitted
    }
    rouder_e = get_values("Rouder half-Cauchy", e_col)
    rouder_v = get_values("Rouder half-Cauchy", v_col)
    flat = get_values(FLAT, e_col)

    same_rouder = np.allclose(
        rouder_e, rouder_v, rtol=0, atol=1e-12
    )
    all_values = np.concatenate(
        [a for pair in fitted_values.values() for a in pair]
        + [rouder_e, rouder_v, flat]
    )
    pad = .08 * max(float(np.ptp(all_values)), 1.)
    limits = (
        all_values.min() - pad,
        all_values.max() + pad,
    )

    fig, axes = plt.subplots(
        2, 3, figsize=(7.8, 8.3), sharex=True, sharey=True
    )
    fig.subplots_adjust(
        left=.255, right=.98, top=.91, bottom=.115,
        wspace=.16, hspace=.22,
    )

    for i, (prior, ax) in enumerate(zip(fitted, axes.flat)):
        e, v = fitted_values[prior]

        # Alternating backgrounds make lineage rows easier to follow.
        for j in range(0, len(names), 2):
            ax.axhspan(
                j-.5, j+.5, color="#F5F7F8", zorder=0
            )

        # Connect only the two fitted-prior methods within each row.
        ax.hlines(
            y, e, v, color="#AAB3B9", lw=1.1, zorder=1
        )

        # Hollow comparator markers remain visible near other points.
        ax.scatter(
            flat, y, marker="^", s=40, facecolors="none",
            edgecolors="#66747D", linewidths=.9, zorder=2,
        )
        ax.scatter(
            rouder_e, y, marker="s", s=23, facecolors="none",
            edgecolors="#C66A26", linewidths=1., zorder=3,
        )
        if not same_rouder:
            ax.scatter(
                rouder_v, y, marker="+", s=25,
                color="#9B5B32", linewidths=.9, zorder=4,
            )

        ax.scatter(
            e, y, marker="o", s=16, color=EVILLE, zorder=5
        )
        ax.scatter(
            v, y, marker="x", s=21, color=VILLE,
            linewidths=.9, zorder=6,
        )

        ax.set_xlim(*limits)
        ax.set_ylim(len(names)-.5, -.5)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.tick_params(
            axis="y", length=0, labelleft=(i % 3 == 0)
        )
        ax.tick_params(axis="x", labelbottom=True)
        ax.xaxis.set_major_locator(
            MaxNLocator(4, prune="both")
        )
        ax.grid(
            axis="x", color="#E4E8EB", lw=.5, zorder=0
        )
        ax.spines["left"].set_visible(False)
        ax.set_title(
            f"({chr(97+i)}) {PRIOR_LABELS[prior]}",
            loc="left", pad=7,
        )

    handles = [
        Line2D(
            [], [], ls="", marker="o", color=EVILLE, ms=4,
            label="Fitted prior: eVille",
        ),
        Line2D(
            [], [], ls="", marker="x", color=VILLE, ms=4,
            label="Fitted prior: Ville",
        ),
        Line2D(
            [], [], ls="", marker="^", color="#66747D",
            markerfacecolor="none", ms=6,
            label="Lebesgue: eVille",
        ),
        Line2D(
            [], [], ls="", marker="s", color="#C66A26",
            markerfacecolor="none", ms=4.5,
            label=(
                "Rouder eVille"
                if same_rouder else "Rouder: eVille"
            ),
        ),
    ]
    if not same_rouder:
        handles.append(
            Line2D(
                [], [], ls="", marker="+", color="#9B5B32",
                ms=5, label="Rouder: Ville",
            )
        )

    # Use the sixth panel for the common legend.
    key = axes.flat[-1]
    key.axis("off")
    key.legend(
        handles=handles, loc="center", frameon=False,
        fontsize=7.5, labelspacing=1.3, handletextpad=.7,
    )

    fig.suptitle(
        "Stopping times across held-out lineages",
        x=.025, y=.98, ha="left", fontsize=11,
    )
    fig.text(
        .025, .948,
        f"{settings['gene']}; recorded Sex = {settings['sex']}; "
        f"{settings['permutations']} paired random orders per lineage",
        fontsize=8,
    )
    fig.text(
        .62, .065,
        "Mean first rejection time (observations)",
        ha="center", fontsize=8,
    )
    fig.text(
        .025, .016,
        "Lower is faster. Numbers in parentheses are target "
        "sample sizes; only N >= 10 is shown.\n"
        "Segments join fitted-prior eVille and Ville for the same "
        "lineage. Horseshoe uses marginal maximum likelihood.",
        fontsize=6.8, linespacing=1.5,
    )

    return _save(
        fig, out, "raw_stopping_times_all_lineages.pdf"
    )


def export_paper_stopping_table(summary, folder, example_target):
    fitted = [
        "Gaussian", "Cauchy", "Laplace",
        "Student-t5", "Horseshoe",
    ]
    benchmarks = ["Improper flat", "Rouder half-Cauchy"]

    e_col = "eville_mean_first_rejection_time_if_all_orders_reject"
    v_col = "ville_mean_first_rejection_time_if_all_orders_reject"

    # Exclude lineages with fewer than ten observations.
    main = summary.loc[summary["target_n"].ge(10)].copy()

    panels = [
        ("All eligible lineages", main),
        (example_target, main.loc[main["target_lineage"].eq(example_target)]),
    ]

    rows = []
    for panel, data in panels:
        if data.empty:
            raise ValueError(f"No eligible results for {panel}.")

        e = data.pivot(
            index="target_lineage", columns="prior", values=e_col
        ).reindex(columns=fitted + benchmarks)

        v = data.pivot(
            index="target_lineage", columns="prior", values=v_col
        ).reindex(index=e.index, columns=fitted)

        # These means require rejection in every random order.
        if not (
            np.isfinite(e.to_numpy()).all()
            and np.isfinite(v.to_numpy()).all()
        ):
            raise ValueError(
                f"{panel}: missing results or at least one nonrejection. "
                "Cannot report these as uncensored mean rejection times."
            )

        # Each lineage receives equal weight.
        for prior in fitted:
            rows.append({
                "Panel": panel,
                "Fitted prior": prior,
                "Fitted prior: eVille": e[prior].mean(),
                "Fitted prior: Ville": v[prior].mean(),
                "Lebesgue: eVille": e["Improper flat"].mean(),
                "Rouder eVille": e["Rouder half-Cauchy"].mean(),
            })

    destination = Path(folder) / "paper_stopping_times.csv"
    pd.DataFrame(rows).to_csv(
        destination, index=False, float_format="%.2f"
    )
    print(f"Saved: {destination}")


def make_figures(results_dir) -> list[Path]:
    """Read consolidated results and create the six PDF figures."""
    folder = Path(results_dir)
    settings = json.loads((folder / "settings.json").read_text())
    summary = pd.read_csv(folder / "summary.csv")
    export_paper_stopping_table(summary, folder, settings["example_target"])
    cohort = pd.read_csv(folder / "cohort.csv")
    paths = pd.read_csv(folder / "example_paths.csv")
    fits = pd.read_csv(folder / "prior_fits.csv")
    if summary.duplicated(["target_lineage", "prior"]).any():
        raise ValueError("summary.csv must have one row per target/prior")
    if set(summary.prior) != set(PRIORS):
        raise ValueError("Expected all seven priors in summary.csv")
    out = folder / "figures"
    out.mkdir(parents=True, exist_ok=True)
    configure_style()
    return [raw_stopping_times_all_lineages(summary, settings, out),
            paired_savings_all_lineages(summary, settings, out),
            conflict_and_flat_baseline(summary, settings, out, fits),
            conflict_and_flat_baseline(
                summary, settings, out, fits,
                lineage=settings["example_target"],
                filename="conflict_and_flat_baseline_example.pdf",
    ),
        example_forest_and_evidence(cohort, paths, fits, settings, out),
        example_running_pvalues(paths, settings, out)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path,
                        default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()
    for path in make_figures(args.results):
        print(path)


if __name__ == "__main__":
    main()

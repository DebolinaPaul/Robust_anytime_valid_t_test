"""Data preparation and fixed-sample descriptive summaries.

Read the chosen Chronos gene column, filter recorded sex, and retain one
ModelID per PatientID before constructing disjoint source/target lineages.
Sources contribute equally weighted standardized lineage means. Confidence
intervals are marginal fixed-sample Gaussian-model intervals, not confidence
sequences and not simultaneous intervals across lineages.
"""
from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd

from effect_intervals import standardized_mean_ci


def load_gene_data(gene_effect_path, models_path, gene, sex="Male"):
    """Read the gene, restrict recorded Sex, then deduplicate patients.

    Among matching rows with a finite score and complete patient/lineage metadata, retain
    the lexicographically first ModelID for each PatientID. This deterministic
    rule is applied globally, before assigning source and target lineages.
    """
    if sex not in ("Male", "Female"):
        raise ValueError(
            "sex must be Male or Female, matching the recorded Sex field")
    header = list(pd.read_csv(gene_effect_path, nrows=0).columns)
    hits = [c for c in header if c == gene or re.match(
        r"^" + re.escape(gene) + r"\s*\(", c)]
    if len(hits) != 1:
        raise ValueError(
            f"Expected exactly one column for {gene!r}; found {hits}.")
    score_column = hits[0]
    id_column = "ModelID" if "ModelID" in header else header[0]
    effects = pd.read_csv(gene_effect_path, usecols=[id_column, score_column])
    effects = effects.rename(
        columns={id_column: "ModelID", score_column: "Chronos"})
    effects["Chronos"] = pd.to_numeric(effects["Chronos"], errors="coerce")
    models = pd.read_csv(models_path, low_memory=False)
    required = ["ModelID", "PatientID", "OncotreeLineage", "Sex"]
    missing = [c for c in required if c not in models]
    if missing:
        raise ValueError(f"Model.csv is missing required columns: {missing}")
    if models["ModelID"].duplicated().any():
        raise ValueError("Model.csv contains duplicate ModelID values.")
    if effects["ModelID"].duplicated().any():
        raise ValueError(
            "The gene-effect file contains duplicate ModelID values.")
    extra = [c for c in ["CellLineName", "OncotreePrimaryDisease"] if c in models]
    data = effects.merge(models[required + extra],
                         on="ModelID", how="left", validate="one_to_one")
    audit = {"gene": gene, "gene_column": score_column,
             "gene_effect_rows": len(effects)}
    normalized_sex = data["Sex"].fillna(
        "").astype(str).str.strip().str.casefold()
    matches_sex = normalized_sex.eq(sex.casefold())
    audit.update(sex_filter=sex, sex_metadata_column="Sex",
                 sex_counts_before_filter={str(k): int(v) for k, v in
                                           data.Sex.fillna("Missing").value_counts().items()},
                 rows_excluded_by_sex_filter=int((~matches_sex).sum()),
                 excluded_sex_counts=data.loc[~matches_sex, "Sex"].fillna(
                     "Missing").value_counts().to_dict(),
                 unknown_or_other_sex_rows_excluded=int((~normalized_sex.isin(["male", "female"])).sum()))
    finite = np.isfinite(data["Chronos"])
    patient_ok = data["PatientID"].notna() & data["PatientID"].fillna(
        "").astype(str).str.strip().ne("")
    lineage_ok = data["OncotreeLineage"].notna() & data["OncotreeLineage"].fillna(
        "").astype(str).str.strip().ne("")
    cancer = ~data["OncotreeLineage"].fillna("").astype(str).str.replace(
        r"[^A-Za-z]", "", regex=True).str.lower().eq("noncancerous")
    audit.update(nonfinite_score_rows=int((~finite).sum()), missing_patient_rows=int((~patient_ok).sum()),
                 missing_lineage_rows=int((~lineage_ok).sum()), noncancerous_rows=int((~cancer).sum()))
    audit["otherwise_eligible_rows_excluded_by_sex_filter"] = int(
        (finite & patient_ok & lineage_ok & cancer & ~matches_sex).sum())
    data = data.loc[finite & patient_ok &
                    lineage_ok & cancer & matches_sex].copy()
    data["Sex"] = sex
    data["ModelID"] = data["ModelID"].astype(str)
    data["PatientID"] = data["PatientID"].astype(str)
    data = data.sort_values("ModelID", kind="stable")
    audit["duplicate_patient_rows_removed"] = int(
        data["PatientID"].duplicated().sum())
    data = data.drop_duplicates(
        "PatientID", keep="first").reset_index(drop=True)
    audit["eligible_distinct_patients"] = len(data)
    if not data.Sex.eq(sex).all():
        raise AssertionError(
            "Every model must match the requested recorded Sex.")
    return data, audit


def summarize(values):
    y = np.asarray(values, float)
    n = len(y)
    mean = float(np.mean(y))
    sd = float(np.std(y, ddof=1)) if n > 1 else float("nan")
    theta = mean / sd if sd > 0 else float("nan")
    theta_se = math.sqrt(1 / n + theta * theta / (2 * (n - 1))
                         ) if n > 1 and sd > 0 else float("nan")
    lower, upper = standardized_mean_ci(
        math.sqrt(n)*theta, n) if np.isfinite(theta) else (float('nan'), float('nan'))
    return {"n": n, "mean_chronos": mean, "sd_chronos": sd,
            "theta_hat": theta, "theta_se_approx": theta_se,
            "theta_ci95_low_approx": theta - 1.96 * theta_se,
            "theta_ci95_high_approx": theta + 1.96 * theta_se,
            "theta_ci95_low": lower, "theta_ci95_high": upper,
            "interval_method": "noncentral-t inversion; fixed-sample Gaussian model"}


def partition(data, target, min_source_n):
    available = sorted(data["OncotreeLineage"].unique())
    if target not in available:
        raise ValueError(
            f"Target lineage {target!r} not present. Available: {available}")
    target_data = data.loc[data.OncotreeLineage.eq(
        target)].sort_values("ModelID").copy()
    if len(target_data) < 2:
        raise ValueError(
            "The target lineage must contain at least two distinct patients for n0=2.")
    source_summaries = []
    excluded = []
    eligible = []
    for lineage, group in data.loc[~data.OncotreeLineage.eq(target)].groupby("OncotreeLineage", sort=True):
        record = {"lineage": lineage, **summarize(group.Chronos)}
        if len(group) < min_source_n:
            excluded.append(
                {**record, "reason": f"fewer than {min_source_n} distinct patients"})
        elif not np.isfinite(record["theta_hat"]):
            excluded.append(
                {**record, "reason": "zero or undefined sample standard deviation"})
        else:
            source_summaries.append(record)
            eligible.append(lineage)
    if len(source_summaries) < 4:
        raise ValueError(
            "At least four eligible source lineages are needed to fit the priors.")
    sources = data.loc[data.OncotreeLineage.isin(eligible)].sort_values(
        ["OncotreeLineage", "ModelID"]).copy()
    if set(sources.PatientID) & set(target_data.PatientID):
        raise AssertionError(
            "Source and target PatientID sets must be disjoint.")
    return target_data, sources, pd.DataFrame(source_summaries), pd.DataFrame(excluded)


def assess_cohort(data, requested_targets, min_source_n=10, min_target_n=10):
    """Account for every retained lineage before fitting or target inspection.

    Eligibility depends on sample size and whether the sample standard
    deviation exists, never on a rejection time or estimated effect magnitude.
    """
    requested = set(requested_targets)
    records = []
    for lineage, group in data.groupby("OncotreeLineage", sort=True):
        count = len(group)
        sd = float(group.Chronos.std(ddof=1)) if count >= 2 else math.nan
        sd_ok = np.isfinite(sd) and sd > 0
        source_eligible = count >= min_source_n and sd_ok
        if lineage not in requested:
            status, reason = "not_requested", "Not included in the explicit --targets selection"
        elif count < 2:
            status, reason = "not_testable", "n0=2 requires at least two target observations"
        elif count < min_target_n:
            status, reason = "excluded_by_target_minimum", f"fewer than {min_target_n} retained target observations"
        elif not sd_ok:
            status, reason = "not_testable", "zero or undefined target sample standard deviation"
        else:
            status, reason = "eligible", ""
        records.append(dict(target_lineage=lineage, target_n=count,
                            source_eligible=bool(source_eligible),
                            requested_as_target=lineage in requested,
                            target_eligible=status == "eligible", status=status,
                            reason=reason, min_source_n=min_source_n,
                            min_target_n=min_target_n))
    return pd.DataFrame(records)

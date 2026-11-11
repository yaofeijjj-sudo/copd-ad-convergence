"""Donor-level association of SEA-AD Axis A UCell scores with DFC CPS_Local."""

from __future__ import annotations

import argparse
import hashlib
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def zscore(series: pd.Series) -> pd.Series:
    sd = series.std(ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError(f"Cannot standardize {series.name}: sd={sd}")
    return (series - series.mean()) / sd


def fit_one(part: pd.DataFrame, supertype: str) -> tuple[dict, pd.DataFrame]:
    frame = part.copy().sort_values("Donor ID").reset_index(drop=True)
    for source, target in [
        ("cps_local", "CPS_Local_z"),
        ("age", "age_z"),
        ("pmi", "PMI_z"),
        ("rin", "RIN_z"),
    ]:
        frame[target] = zscore(frame[source].astype(float))

    if set(frame["sex"].unique()) - {"Female", "Male"}:
        raise ValueError(f"Unexpected sex labels for {supertype}: {sorted(frame['sex'].unique())}")
    if set(frame["study"].unique()) - {"ACT", "ADRC Clinical Core"}:
        raise ValueError(f"Unexpected study labels for {supertype}: {sorted(frame['study'].unique())}")

    frame["sex_Male"] = frame["sex"].eq("Male").astype(float)
    frame["study_ADRC_Clinical_Core"] = frame["study"].eq("ADRC Clinical Core").astype(float)
    xcols = [
        "CPS_Local_z",
        "age_z",
        "sex_Male",
        "PMI_z",
        "RIN_z",
        "study_ADRC_Clinical_Core",
    ]
    x = sm.add_constant(frame[xcols], has_constant="add").astype(float)
    y = frame["AxisA_full_measured_UCell"].astype(float)
    rank = int(np.linalg.matrix_rank(x.to_numpy()))
    if rank != x.shape[1]:
        raise ValueError(f"Design is rank deficient for {supertype}: rank={rank}, columns={x.shape[1]}")

    fitted = sm.OLS(y, x).fit(cov_type="HC3")
    beta = float(fitted.params["CPS_Local_z"])
    se = float(fitted.bse["CPS_Local_z"])
    pvalue = float(fitted.pvalues["CPS_Local_z"])
    ci_low, ci_high = map(float, fitted.conf_int().loc["CPS_Local_z"])
    row = {
        "Supertype": supertype,
        "n_donors": int(frame["Donor ID"].nunique()),
        "n_nuclei_total": int(frame["n_nuclei"].sum()),
        "median_nuclei_per_donor": float(frame["n_nuclei"].median()),
        "median_donor_UCell": float(y.median()),
        "CPS_Local_mean": float(frame["cps_local"].mean()),
        "CPS_Local_sd": float(frame["cps_local"].std(ddof=1)),
        "beta_UCell_per_1SD_CPS_Local": beta,
        "HC3_se": se,
        "ci95_low": ci_low,
        "ci95_high": ci_high,
        "p_two_sided": pvalue,
        "model_r_squared": float(fitted.rsquared),
        "design_rows": int(x.shape[0]),
        "design_columns": int(x.shape[1]),
        "design_rank": rank,
        "primary_endpoint": supertype == "Micro-PVM_2",
        "model": "OLS with HC3 SE: UCell ~ CPS_Local_z + age_z + sex + PMI_z + RIN_z + study",
    }
    keep = [
        "Donor ID",
        "Supertype",
        "AxisA_full_measured_UCell",
        "n_nuclei",
        "cps_local",
        "CPS_Local_z",
        "age",
        "sex",
        "pmi",
        "rin",
        "study",
    ]
    return row, frame[keep]

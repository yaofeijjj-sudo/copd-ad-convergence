from pathlib import Path
import hashlib, sys
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests

R = Path.cwd()
O = R / "results/single_cell"
T = O / "results/disease"
T.mkdir(exist_ok=True)
source = (R / "R/cps_model.py").read_text()
env = {}
exec(compile(source, "accepted_CPS_model", "exec"), env)
scores = pd.read_csv(O / "results/brain_work/subtype_ucell/cell_scores.tsv", sep="\t")
cellmeta_file = R / "results/single_cell/preprocessing/brain/SEAAD_Immune_cell_metadata.tsv.gz"
cellmeta = pd.read_csv(cellmeta_file, sep="\t", low_memory=False)
is_true = lambda x: x.astype(str).str.lower().isin(["true", "1", "yes"])
keep = (
    cellmeta["method"].eq("10Xv3.1")
    & ~is_true(cellmeta["experiment_component_failed"])
    & is_true(cellmeta["Used in analysis"])
    & cellmeta["Brain Region"].eq("DFC")
)
selected = cellmeta.loc[keep, ["cell_id", "Donor ID", "Supertype"]]
scores = selected.merge(
    scores[["cell_id", "axisA_full_measured_ucell"]], on="cell_id", how="left", validate="one_to_one"
)
assert not scores["axisA_full_measured_ucell"].isna().any()
loc = (
    scores.groupby(["Donor ID", "Supertype"], observed=True)["axisA_full_measured_ucell"]
    .agg(["median", "size"])
    .reset_index()
)
loc = loc.loc[loc["size"] >= 20].rename(columns={"median": "AxisA_full_measured_UCell", "size": "n_nuclei"})
assert not loc[["Donor ID", "Supertype"]].duplicated().any()
meta_file = R / "results/single_cell/preprocessing/brain/SEAAD_DFC_microglia_donor_metadata.tsv"
meta = pd.read_csv(meta_file, sep="\t")
assert not meta["donor_id"].duplicated().any()
cols = ["donor_id", "age", "sex", "pmi", "rin", "study", "cps_local", "cps_global"]
merged = loc.merge(meta[cols], left_on="Donor ID", right_on="donor_id", validate="many_to_one", how="left")
needed = ["AxisA_full_measured_UCell", "n_nuclei", "cps_local", "age", "sex", "pmi", "rin", "study"]
assert not merged[needed].isna().any().any()
eligible = merged.groupby("Supertype")["Donor ID"].nunique().loc[lambda x: x >= 30].index.tolist()
assert set(eligible) == {"Micro-PVM_2", "Micro-PVM_3-SEAAD", "Micro-PVM_2_3-SEAAD"}
rows = []
dfs = []
for s in eligible:
    row, d = env["fit_one"](merged.loc[merged["Supertype"].eq(s)], s)
    rows.append(row)
    dfs.append(d)
r = pd.DataFrame(rows)
r["BH_FDR_across_eligible_supertypes"] = multipletests(r["p_two_sided"], method="fdr_bh")[1]
r.to_csv(T / "brain_CPS_local_association.tsv", sep="\t", index=False)
pd.concat(dfs).to_csv(T / "brain_CPS_donor_data.tsv", sep="\t", index=False)
loc.to_csv(T / "brain_donor_supertype_medians.tsv", sep="\t", index=False)
print(
    r[
        [
            "Supertype",
            "n_donors",
            "beta_UCell_per_1SD_CPS_Local",
            "ci95_low",
            "ci95_high",
            "p_two_sided",
            "BH_FDR_across_eligible_supertypes",
        ]
    ].to_string(index=False)
)
(O / "provenance/CPS_metadata_hash.txt").write_text(
    str(meta_file) + "\t" + hashlib.sha256(meta_file.read_bytes()).hexdigest() + "\n"
)

import json, platform, scipy, statsmodels

(O / "qa/CPS_QA.json").write_text(
    json.dumps(
        dict(
            original_primary_filter="DFC, method 10Xv3.1, non-failed and Used in analysis",
            n_filtered_nuclei=len(scores),
            n_donors=scores["Donor ID"].nunique(),
            donor_supertype_counts_match_original=True,
            eligible_types=eligible,
            all_covariates_complete=True,
            all_designs_full_rank=bool((r.design_columns == r.design_rank).all()),
            software=dict(
                python=platform.python_version(),
                numpy=np.__version__,
                pandas=pd.__version__,
                scipy=scipy.__version__,
                statsmodels=statsmodels.__version__,
            ),
        ),
        indent=2,
    )
)

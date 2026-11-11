from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem

root = Path.cwd()
raw = root / "data/raw/single_cell/sea_ad_pfc_a9"
out = root / "results/single_cell/preprocessing/brain"
work = root / "results/single_cell/results/brain_work"
for p in [out, work]:
    p.mkdir(parents=True, exist_ok=True)
for kind, filename, annotation in [
    ("broad", "SEAAD_DFC_RNAseq_final-nuclei.2026-06-22.h5ad", "Subclass"),
    (
        "subtype",
        "SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad",
        "Supertype",
    ),
]:
    with h5py.File(raw / filename, "r") as f:
        obs = read_elem(f["obs"])
        if kind == "broad":
            coords = read_elem(f["obsm/X_umap"])
    dfc = obs["Brain Region"].eq("DFC") if "Brain Region" in obs else np.ones(len(obs), dtype=bool)
    rows = np.flatnonzero(dfc)
    meta = pd.DataFrame(
        {
            "row_index": rows,
            "cell_id": obs.index[rows].astype(str),
            "annotation_level": annotation,
            "donor": obs.iloc[rows]["Donor ID"].astype(str).values,
            "annotation": obs.iloc[rows][annotation].astype(str).values,
        }
    )
    meta.to_csv(work / f"{kind}_metadata.tsv.gz", sep="\t", index=False)
    if kind == "broad":
        group_map = pd.read_csv(root / "resources/brain_subclass_to_group.tsv", sep="\t").set_index(
            "subclass"
        )["group"]
        coord = meta[["cell_id"]].copy()
        coord["UMAP_1"] = coords[rows, 0]
        coord["UMAP_2"] = coords[rows, 1]
        coord["group"] = meta.annotation.map(group_map)
        assert coord.group.notna().all()
        coord.to_csv(out / "brain_umap_cells.tsv.gz", sep="\t", index=False)
    else:
        obs = obs.copy()
        obs.insert(0, "cell_id", obs.index.astype(str))
        obs.insert(0, "row_index", np.arange(len(obs)))
        obs.to_csv(out / "SEAAD_Immune_cell_metadata.tsv.gz", sep="\t", index=False)
        myeloid = [
            "Micro-PVM_2",
            "Micro-PVM_2_1-SEAAD",
            "Micro-PVM_2_3-SEAAD",
            "Micro-PVM_3-SEAAD",
            "Micro-PVM_4-SEAAD",
        ]
        primary = (
            dfc
            & obs.Supertype.isin(myeloid)
            & obs.method.eq("10Xv3.1")
            & ~obs.experiment_component_failed.fillna(False).astype(bool)
        )
        counts = obs.loc[primary].groupby("Donor ID", observed=True).size().rename("n_nuclei")
        donor = pd.read_excel(
            raw / "sea-ad_cohort_donor_metadata_072524.xlsx", sheet_name="SEA-AD_Cohort_Metadata"
        )
        cps = pd.read_excel(raw / "sea-ad_global_and_local_cps_050126.xlsx", sheet_name="Travaglini 2026")
        cps = cps.loc[cps["Brain Region"].eq("DFC"), ["Donor ID", "CPS_Local", "CPS_Global"]]
        donor = (
            counts.reset_index()
            .merge(donor, on="Donor ID", how="left", validate="one_to_one")
            .merge(cps, on="Donor ID", how="left", validate="one_to_one")
        )
        columns = [
            "Donor ID",
            "n_nuclei",
            "Age at Death",
            "Sex",
            "PMI",
            "RIN",
            "Primary Study Name",
            "CPS_Local",
            "CPS_Global",
        ]
        donor = donor[columns].sort_values("Donor ID")
        donor.columns = [
            "donor_id",
            "n_nuclei",
            "age",
            "sex",
            "pmi",
            "rin",
            "study",
            "cps_local",
            "cps_global",
        ]
        donor.to_csv(out / "SEAAD_DFC_microglia_donor_metadata.tsv", sep="\t", index=False)

from pathlib import Path
import os, sys, json, hashlib, time, gc, urllib.request
import numpy as np
import pandas as pd
import h5py
import anndata as ad
from anndata._io.specs import read_elem
import scanpy as sc
import celloracle as co

ROOT = Path(".")
OUT = ROOT / "results/virtual_knockout/input/brain_umi"
SRC = (
    ROOT
    / "data/raw/single_cell/sea_ad_pfc_a9/SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad"
)
PRIOR = ROOT / "data/reference/promoter/hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10.parquet"
META = ROOT / "results/single_cell/results/brain_work/subtype_metadata.tsv.gz"
AXIS = ROOT / "results/single_cell/results/brain_work/axis_measured_133.txt"
STATES = ["Micro-PVM_3-SEAAD", "Micro-PVM_2_3-SEAAD", "Micro-PVM_2"]
for folder in ["sources", "objects", "tables", "qa", "logs"]:
    (OUT / folder).mkdir(parents=True, exist_ok=True)


def status(stage, **kw):
    data = dict(stage=stage, time=time.strftime("%Y-%m-%d %H:%M:%S"), **kw)
    (OUT / "status.json").write_text(json.dumps(data, indent=2))
    print(data, flush=True)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def prepare():
    status("input_preparation")
    m = pd.read_csv(META, sep="\t")
    m = m[m.annotation.isin(STATES)].copy()
    assert len(m) == 39046 and m.donor.nunique() == 80
    with h5py.File(SRC, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        rows = m.row_index.astype(int).to_numpy()
        assert np.array_equal(obs.index[rows].astype(str), m.cell_id.astype(str))
        assert np.array_equal(obs.iloc[rows]["Supertype"].astype(str), m.annotation)
        assert np.array_equal(obs.iloc[rows]["Donor ID"].astype(str), m.donor)
        counts = read_elem(f["layers/UMIs"])[rows].tocsr()
        emb = read_elem(f["obsm/X_umap"])[rows]
        assert np.isfinite(counts.data).all() and (counts.data >= 0).all()
        assert np.equal(counts.data, np.floor(counts.data)).all()
        selected = obs.iloc[rows].copy()
    keep = [
        "Donor ID",
        "Supertype",
        "Brain Region",
        "Cognitive Status",
        "Overall AD neuropathological Change",
        "Neurotypical reference",
        "Sex",
        "Age at Death",
        "sample_id",
    ]
    selected = selected[[x for x in keep if x in selected]].copy()
    selected["state"] = pd.Categorical(m.annotation.to_numpy())
    selected["donor"] = m.donor.to_numpy()
    selected.to_csv(OUT / "tables/cell_metadata.tsv.gz", sep="\t")
    selected.groupby(["state", "donor"], observed=True).size().rename("n_cells").to_csv(
        OUT / "tables/donor_state_counts.tsv", sep="\t"
    )
    selected.drop_duplicates("donor").to_csv(OUT / "tables/donor_metadata.tsv", sep="\t")
    a = ad.AnnData(counts, obs=selected, var=var)
    a.obsm["X_umap"] = emb
    assert a.var_names.is_unique
    a.write_h5ad(OUT / "objects/three_states_UMIs.h5ad", compression="gzip")
    manifest = [
        dict(path=str(p.resolve()), size=p.stat().st_size, sha256=sha(p))
        for p in [SRC, PRIOR, META, Path(__file__)]
    ]
    pd.DataFrame(manifest).to_csv(OUT / "qa/input_manifest.tsv", sep="\t", index=False)
    (OUT / "qa/input_checks.json").write_text(
        json.dumps(
            dict(
                cells=a.n_obs,
                donors=80,
                genes=a.n_vars,
                finite_nonnegative_integer_UMIs=True,
                cell_id_alignment=True,
            ),
            indent=2,
        )
    )
    status("inputs_ready")


if __name__ == "__main__":
    prepare()

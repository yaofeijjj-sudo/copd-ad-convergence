"""Annotate QC-passing GSE302339 cells with the CellTypist Human_Lung_Atlas model."""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import anndata
import celltypist
import numpy as np
import pandas as pd
import scanpy as sc


PROJECT_ROOT = Path(".")
RESULT_DIR = PROJECT_ROOT / "results/single_cell/preprocessing" / "lung"
H5_DIR = PROJECT_ROOT / "data/raw/lung_h5"
MODEL_FILE = RESULT_DIR / "inputs/Human_Lung_Atlas.pkl"
QC_FILE = RESULT_DIR / "objects/cell_qc_all.tsv.gz"
OUT_DIR = RESULT_DIR / "objects/celltypist_by_library"
OUT_DIR.mkdir(parents=True, exist_ok=True)


MYELOID = {
    "Alveolar Mph CCL3+",
    "Alveolar Mph MT-positive",
    "Alveolar Mph proliferating",
    "Alveolar macrophages",
    "Classical monocytes",
    "DC1",
    "DC2",
    "Interstitial Mph perivascular",
    "Mast cells",
    "Migratory DCs",
    "Monocyte-derived Mph",
    "Non-classical monocytes",
    "Plasmacytoid DCs",
}
T_NK = {"CD4 T cells", "CD8 T cells", "NK cells", "T cells proliferating"}
B_PLASMA = {"B cells", "Plasma cells"}
ENDOTHELIAL = {
    "EC aerocyte capillary",
    "EC arterial",
    "EC general capillary",
    "EC venous pulmonary",
    "EC venous systemic",
    "Lymphatic EC differentiating",
    "Lymphatic EC mature",
    "Lymphatic EC proliferating",
}
STROMAL = {
    "Adventitial fibroblasts",
    "Alveolar fibroblasts",
    "Mesothelium",
    "Myofibroblasts",
    "Peribronchial fibroblasts",
    "Pericytes",
    "SM activated stress response",
    "Smooth muscle",
    "Smooth muscle FAM83D+",
    "Subpleural fibroblasts",
}
OTHER = {"Hematopoietic stem cells"}


def broad_label(fine: str) -> str:
    if fine in MYELOID:
        return "Myeloid"
    if fine in T_NK:
        return "T/NK"
    if fine in B_PLASMA:
        return "B/Plasma"
    if fine in ENDOTHELIAL:
        return "Endothelial"
    if fine in STROMAL:
        return "Stromal"
    if fine in OTHER:
        return "Other/Progenitor"
    return "Epithelial"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


qc = pd.read_csv(
    QC_FILE,
    sep="\t",
    usecols=["barcode", "gsm", "passes_QC"],
    dtype={"barcode": str, "gsm": str, "passes_QC": bool},
)
qc = qc.loc[qc["passes_QC"]].copy()
model = celltypist.models.Model.load(model=str(MODEL_FILE))

summaries: list[dict] = []
for idx, gsm in enumerate(sorted(qc["gsm"].unique()), start=1):
    out_file = OUT_DIR / f"{gsm}_celltypist.tsv.gz"
    if out_file.exists():
        ann = pd.read_csv(out_file, sep="\t")
        print(f"[{idx}/65] {gsm} checkpoint found", flush=True)
    else:
        h5_candidates = sorted(H5_DIR.glob(f"{gsm}_*.h5"))
        if len(h5_candidates) != 1:
            raise RuntimeError(f"Expected one H5 file for {gsm}, found {len(h5_candidates)}")
        print(f"[{idx}/65] {gsm} reading and annotating", flush=True)
        adata = sc.read_10x_h5(h5_candidates[0])
        adata.var_names_make_unique()
        keep = qc.loc[qc["gsm"] == gsm, "barcode"].to_numpy()
        missing = np.setdiff1d(keep, adata.obs_names.to_numpy())
        if missing.size:
            raise RuntimeError(f"{gsm}: {missing.size} QC barcodes absent from H5")
        adata = adata[keep].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        prediction = celltypist.annotate(
            adata,
            model=model,
            mode="best match",
            majority_voting=False,
        )
        labels = prediction.predicted_labels["predicted_labels"].astype(str)
        probabilities = prediction.probability_matrix.max(axis=1)
        ann = pd.DataFrame(
            {
                "barcode": labels.index.astype(str),
                "gsm": gsm,
                "celltypist_label": labels.to_numpy(),
                "celltypist_conf_score": probabilities.to_numpy(),
            }
        )
        ann["broad_celltype"] = ann["celltypist_label"].map(broad_label)
        ann.to_csv(out_file, sep="\t", index=False, compression="gzip")
    summaries.append(
        {
            "gsm": gsm,
            "n_annotated": len(ann),
            "n_fine_labels": ann["celltypist_label"].nunique(),
            "median_conf_score": ann["celltypist_conf_score"].median(),
        }
    )

summary = pd.DataFrame(summaries)
summary.to_csv(RESULT_DIR / "tables/celltypist_library_summary.tsv", sep="\t", index=False)

all_annotations = pd.concat(
    [pd.read_csv(p, sep="\t") for p in sorted(OUT_DIR.glob("*_celltypist.tsv.gz"))],
    ignore_index=True,
)
all_annotations.to_csv(
    RESULT_DIR / "objects/celltypist_all.tsv.gz", sep="\t", index=False, compression="gzip"
)
counts = (
    all_annotations.groupby(["broad_celltype", "celltypist_label"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
    .sort_values(["broad_celltype", "n_cells"], ascending=[True, False])
)
counts.to_csv(RESULT_DIR / "tables/celltypist_label_counts.tsv", sep="\t", index=False)

run_info = {
    "python": platform.python_version(),
    "celltypist": celltypist.__version__,
    "scanpy": sc.__version__,
    "anndata": anndata.__version__,
    "model": MODEL_FILE.name,
    "model_sha256": sha256(MODEL_FILE),
    "mode": "best match",
    "majority_voting": False,
    "normalization": "scanpy.pp.normalize_total(target_sum=1e4), then scanpy.pp.log1p",
    "n_cells": int(len(all_annotations)),
    "n_fine_labels": int(all_annotations["celltypist_label"].nunique()),
}
(RESULT_DIR / "qa/celltypist_run_info.json").write_text(json.dumps(run_info, indent=2) + "\n")
print(json.dumps(run_info, indent=2), flush=True)

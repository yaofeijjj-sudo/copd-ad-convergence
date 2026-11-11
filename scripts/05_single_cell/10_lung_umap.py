from pathlib import Path
import sys

sys.path.insert(0, str(Path.cwd() / "R"))
from runtime import params
import os

os.environ["NUMBA_NUM_THREADS"] = "6"
import scanpy as sc, anndata as ad, pandas as pd, numpy as np, json, hashlib, scipy.sparse as sp, harmonypy, importlib.metadata as md

root = Path(".")
out = root / "results/single_cell/preprocessing_figures/lung_umap_reference"
cache = root / "results/cache/lung_umap"
cache.mkdir(parents=True, exist_ok=True)
for name in ["source_data", "qa"]:
    (out / name).mkdir(parents=True, exist_ok=True)
src = root / "results/single_cell/results/lung/objects/axisA_cell_scores.tsv.gz"
d = pd.read_csv(src, sep="\t")
assert d.cell_id.is_unique and d.passes_QC.all()
obs = d.set_index("cell_id")
mats = []
ids = []
sources = []
genes = None
for i, (gsm, g) in enumerate(d.groupby("gsm", sort=True)):
    f = list((root / "data/raw/lung_h5").glob(gsm + "_*"))[0]
    a = sc.read_10x_h5(f)
    a.var_names_make_unique()
    if genes is None:
        genes = a.var_names.copy()
    assert genes.equals(a.var_names)
    a = a[g.barcode.to_numpy()].copy()
    mats.append(a.X)
    ids.extend(g.cell_id)
    sources.append({"path": str(f), "sha256": hashlib.sha256(f.read_bytes()).hexdigest()})
    if i % 10 == 0:
        print("Loaded libraries", i + 1, flush=True)
a = ad.AnnData(sp.vstack(mats, format="csr"), obs=obs.loc[ids].copy(), var=pd.DataFrame(index=genes))
del mats
print("Matrix", a.shape, flush=True)
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)
sc.pp.highly_variable_genes(a, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key="gsm")
a.var.to_csv(out / "source_data/gene_selection.tsv", sep="\t")
a = a[:, a.var.highly_variable].copy()
print("HVG", a.shape, flush=True)
sc.pp.scale(a, max_value=10)
sc.pp.pca(a, n_comps=50, random_state=params["umap_seed"])
print("PCA complete", flush=True)
ho = harmonypy.run_harmony(
    a.obsm["X_pca"], a.obs, "gsm", max_iter_harmony=50, random_state=params["umap_seed"]
)
z = ho.Z_corr
if z.shape != a.obsm["X_pca"].shape:
    z = z.T
assert z.shape == a.obsm["X_pca"].shape
a.obsm["X_pca_harmony"] = z
np.savez_compressed(cache / "pca_harmony.npz", cell_id=a.obs_names.to_numpy(dtype=str), pca=z)
print("Harmony complete", flush=True)
b = ad.AnnData(sp.csr_matrix((a.n_obs, 1)), obs=a.obs.copy())
b.obsm["X_pca_harmony"] = z
del a, ho
sc.pp.neighbors(b, n_neighbors=15, n_pcs=50, use_rep="X_pca_harmony", random_state=params["umap_seed"])
sc.tl.umap(b, min_dist=0.5, spread=1, random_state=params["umap_seed"])
print("UMAP complete", flush=True)
t = b.obs[["gsm", "barcode", "subject_id", "broad_celltype", "AxisA_UCell"]].copy()
t["UMAP_1"] = b.obsm["X_umap"][:, 0]
t["UMAP_2"] = b.obsm["X_umap"][:, 1]
assert np.isfinite(t[["UMAP_1", "UMAP_2", "AxisA_UCell"]]).all().all()
t.to_csv(out / "source_data/lung_umap_cells.tsv.gz", sep="\t")
b.write_h5ad(cache / "embedding_graph.h5ad", compression="gzip")
qa = {
    "n_cells": len(t),
    "n_donors": t.subject_id.nunique(),
    "n_libraries": t.gsm.nunique(),
    "groups": t.broad_celltype.value_counts().to_dict(),
    "score_range": [float(t.AxisA_UCell.min()), float(t.AxisA_UCell.max())],
    "source_score_sha256": hashlib.sha256(src.read_bytes()).hexdigest(),
    "source_scores": str(src),
    "source_matrices": sources,
    "seed": 0,
    "embedding": True,
    "author_workflow": "scanpy_workflow.zip/1_preprocessing_doublet_detection.ipynb",
    "versions": {x: md.version(x) for x in ["scanpy", "harmonypy", "anndata", "numpy", "umap-learn"]},
}
(out / "qa/data_verification.json").write_text(json.dumps(qa, indent=2))
print("DONE", flush=True)

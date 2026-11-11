# SEA-AD DFC single-nucleus reference for CIBERSORTx: 250 nuclei per broad type, balanced across donors,
# raw UMI counts from layers/UMIs. Writes a Matrix Market file (Ensembl IDs) for 03_brain_reference.R.
# Run from project root:
#   _workspace/venvs/seaad_dfc_axisA_localization_20260831_candidate/bin/python -I 03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018/code/02_brain_reference_extract.py
import h5py, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio, os
H5 = "cell_localization/raw_data/sea_ad_pfc_a9_20260622/SEAAD_DFC_RNAseq_final-nuclei.2026-06-22.h5ad"
OUT = "03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018"
N_PER, SEED = 250, 20261018
f = h5py.File(H5, "r"); obs = f["obs"]
def col(k):
    g = obs[k]
    if isinstance(g, h5py.Group):
        cats = g["categories"][:].astype(str); return cats[g["codes"][:]]
    return g[:].astype(str)
cls, sub, donor = col("Class"), col("Subclass"), col("Donor ID")
broad = np.where(cls == "Neuronal: Glutamatergic", "Excitatory", np.where(cls == "Neuronal: GABAergic", "Inhibitory", sub))
broad = np.array([b.replace(" & ", "_").replace(" ", "_") for b in broad])
keep_types = ["Excitatory", "Inhibitory", "Astrocyte", "Oligodendrocyte", "OPC", "Immune", "Endothelial", "VLMC_Perivascular"]
d = pd.DataFrame({"i": np.arange(len(broad)), "type": broad, "donor": donor})
d = d[d["type"].isin(keep_types)]
rng = np.random.default_rng(SEED)
d["u"] = rng.random(len(d)); d = d.sort_values(["type", "donor", "u"]); d["r"] = d.groupby(["type", "donor"]).cumcount()
d["v"] = rng.random(len(d)); d = d.sort_values(["type", "r", "v"]); sel = d.groupby("type").head(N_PER).sort_values("i")
os.makedirs(f"{OUT}/qa", exist_ok=True); os.makedirs(f"{OUT}/source_data", exist_ok=True)
sel.groupby("type").agg(n=("i", "size"), donors=("donor", "nunique"), max_one_donor=("donor", lambda x: x.value_counts().max())).to_csv(f"{OUT}/qa/brain_reference_donor_balance.tsv", sep="\t")
L = f["layers"]["UMIs"]; indptr = L["indptr"]; idx = sel["i"].to_numpy()
rows_data, rows_ind, ptr = [], [], [0]
for i in idx:
    a, b = indptr[i], indptr[i + 1]
    rows_data.append(L["data"][a:b]); rows_ind.append(L["indices"][a:b]); ptr.append(ptr[-1] + (b - a))
M = sp.csr_matrix((np.concatenate(rows_data), np.concatenate(rows_ind), np.array(ptr)), shape=(len(idx), L.attrs["shape"][1] if "shape" in L.attrs else f["X"].attrs["shape"][1]))
assert np.allclose(M.data, np.round(M.data)), "UMIs layer is not integer counts"
ens = f["var"]["gene_ids"][:].astype(str)
sio.mmwrite(f"{OUT}/source_data/brain_ref_counts.mtx", M.T.tocoo())          # genes x cells
pd.DataFrame({"ensembl": ens}).to_csv(f"{OUT}/source_data/brain_ref_genes.tsv", sep="\t", index=False)
sel[["i", "type", "donor"]].to_csv(f"{OUT}/source_data/brain_ref_cells.tsv", sep="\t", index=False)
print(sel["type"].value_counts().to_dict(), "nuclei", len(sel), "genes", M.shape[1], "nnz", M.nnz)

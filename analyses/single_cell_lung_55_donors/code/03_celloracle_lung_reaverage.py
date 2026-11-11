# CellOracle lung (alveolar macrophages, 30,000 cells): the readout is the per-donor mean of the per-cell delta_X, then the
# equal-weight mean across donors (CellOracle_AxisA_Entrez137_20261001/code/finalize.py). cood8 and copd8 are one donor, so
# the lung donor means are recomputed with cood8 -> copd8 (55 donors). The networks and simulations are unchanged (cells are
# pooled for fitting); only this aggregation changes. Brain values are copied. The ranking follows
# CellOracle_ranking_scale_20261019/code/01_ranking.py unchanged.
# Run from project root: _workspace/envs/celloracle/bin/python -I 03_candidate_analysis/Lung_sc_donor_merge_20261029/code/03_celloracle_lung_reaverage.py
import numpy as np, pandas as pd
C = "03_candidate_analysis/CellOracle_AxisA_Entrez137_20261001/"
O = "03_candidate_analysis/Lung_sc_donor_merge_20261029/results/tables/"
p_all = pd.read_csv(C + "tables/all_common_TF_paired_gene_responses.tsv.gz", sep="\t", low_memory=False, dtype={"entrez_id": str})
base = np.load(C + "objects/lung_KO/baseline_AxisA.npz"); genes = base["genes"].astype(str); donors0 = base["donors"].astype(str)
donors1 = np.where(donors0 == "cood8", "copd8", donors0)

def equal_donor_mean(x, donors):
    d = sorted(set(donors)); return np.array([x[donors == k].mean(axis=0, dtype=np.float64) for k in d]).mean(axis=0), len(d)

b0, n0 = equal_donor_mean(base["expression"], donors0); b1, n1 = equal_donor_mean(base["expression"], donors1)
assert (n0, n1) == (56, 55)
main = p_all[p_all["mode"] == "main"].copy()
rows = []
for tf, g in main.groupby("model_TF_lung", sort=False):
    arr = np.load(C + f"objects/lung_KO/{tf}_main.npz"); assert list(arr["axis_genes"].astype(str)) == list(genes)
    d0, _ = equal_donor_mean(arr["delta_axis"], donors0); d1, _ = equal_donor_mean(arr["delta_axis"], donors1)
    rows.append(pd.DataFrame({"model_TF_lung": tf, "gene_lung": genes, "delta0": d0, "delta1": d1, "base0": b0, "base1": b1}))
new = pd.concat(rows)
m = main.merge(new, on=["model_TF_lung", "gene_lung"], how="left", validate="one_to_one")
inm = m.in_model_lung.astype(str).eq("True")
# gate: the original key reproduces the frozen lung readout
assert m.loc[inm, "delta0"].notna().all()
assert np.allclose(m.loc[inm, "delta0"], m.loc[inm, "lung_delta"], rtol=1e-10, atol=1e-14)
assert np.allclose(m.loc[inm, "base0"], m.loc[inm, "lung_baseline"], rtol=1e-10, atol=1e-14)
m.loc[inm, "lung_delta"] = m.loc[inm, "delta1"]; m.loc[inm, "lung_baseline"] = m.loc[inm, "base1"]
m.loc[inm, "lung_simulated"] = m.loc[inm, "base1"] + m.loc[inm, "delta1"]
self_ = m.forced_self_lung.astype(str).eq("True")
m.loc[inm, "lung_downstream_delta"] = m.loc[inm, "delta1"].where(~self_[inm])
m.drop(columns=["delta0", "delta1", "base0", "base1"]).to_csv(O + "all_common_TF_paired_gene_responses_lung55.tsv.gz", sep="\t", index=False)

# ---- ranking (01_ranking.py) ----
p = m[inm & m.in_model_brain.astype(str).eq("True")]
el = pd.read_csv(C + "tables/all_TF_eligibility_and_execution.tsv.gz", sep="\t", low_memory=False)
deg = el.pivot_table(index="gene", columns="side", values="outgoing_nonzero", aggfunc="first")
out = []
for tf, g in p.groupby("TF"):
    l, b = g.lung_delta.values, g.brain_delta.values
    out.append(dict(TF=tf, n_evaluable=len(g),
        shared_down=int(((l < 0) & (b < 0)).sum()), shared_up=int(((l > 0) & (b > 0)).sum()),
        opposite=int(((l > 0) & (b < 0)).sum() + ((l < 0) & (b > 0)).sum()),
        lung_only_change=int(((l != 0) & (b == 0)).sum()), brain_only_change=int(((l == 0) & (b != 0)).sum()),
        mean_signed_delta_lung=l.mean(), mean_signed_delta_brain=b.mean(),
        mean_abs_delta_lung=np.abs(l).mean(), mean_abs_delta_brain=np.abs(b).mean(),
        degree_lung=deg.loc[tf, "lung"] if tf in deg.index else np.nan, degree_brain=deg.loc[tf, "brain"] if tf in deg.index else np.nan))
r = pd.DataFrame(out)
r["shared_down_frac"] = r.shared_down / r.n_evaluable
r["shared_down_per_100_outdegree_mean"] = 100 * r.shared_down / r[["degree_lung", "degree_brain"]].mean(axis=1)
for s in ["lung", "brain"]:
    r[f"response_percentile_{s}"] = 100 * r[f"mean_abs_delta_{s}"].rank(pct=True)
    r[f"signed_percentile_{s}"] = 100 * r[f"mean_signed_delta_{s}"].rank(pct=True)
r = r.sort_values(["shared_down", "TF"], ascending=[False, True]).reset_index(drop=True)
r.insert(0, "rank_shared_down", r.shared_down.rank(ascending=False, method="min").astype(int))
assert len(r) == 357
r.to_csv(O + "all357_TF_ranking.tsv", sep="\t", index=False)
g = p[p.TF.isin(["STAT3", "CEBPB", "LYL1"])][["TF", "entrez_id", "gene", "lung_baseline", "lung_delta", "brain_baseline", "brain_delta"]].copy()
for s in ["lung", "brain"]:
    g[f"{s}_pct_of_log_baseline_old"] = 100 * g[f"{s}_delta"] / g[f"{s}_baseline"]
    g[f"{s}_linear_pct_change"] = 100 * (np.expm1(g[f"{s}_baseline"] + g[f"{s}_delta"]) - np.expm1(g[f"{s}_baseline"])) / np.expm1(g[f"{s}_baseline"])
g.to_csv(O + "Fig6def_source_signed_logratio.tsv", sep="\t", index=False)

# ---- comparison with the frozen ranking ----
f = pd.read_csv("03_candidate_analysis/CellOracle_ranking_scale_20261019/results/all357_TF_ranking.tsv", sep="\t")
c = f.merge(r, on="TF", suffixes=("_frozen", "_merged"))
keep = ["rank_shared_down", "shared_down", "shared_up", "opposite", "response_percentile_lung"]
show = c[(c.rank_shared_down_frozen <= 12) | (c.rank_shared_down_merged <= 12) | c.TF.isin(["KLF4", "PPARG", "MEF2C", "BHLHE41"])]
print(show[["TF"] + [f"{k}_{s}" for k in keep for s in ("frozen", "merged")]].sort_values("rank_shared_down_merged").to_string(index=False))
lr = lambda x: x.rank(ascending=False, method="min")
for t in ["KLF4", "PPARG", "MEF2C", "BHLHE41"]:
    print(t, "lung rank by mean |delta|: frozen", int(lr(f.mean_abs_delta_lung)[f.TF == t].iloc[0]), "merged", int(lr(r.mean_abs_delta_lung)[r.TF == t].iloc[0]),
          "| brain rank", int(lr(r.mean_abs_delta_brain)[r.TF == t].iloc[0]))
print("changed shared_down:", int((c.shared_down_frozen != c.shared_down_merged).sum()), "of 357; max |change|", int((c.shared_down_frozen - c.shared_down_merged).abs().max()))
print("Spearman shared_down vs mean outdegree:", round(r.shared_down.corr(r[["degree_lung", "degree_brain"]].mean(axis=1), method="spearman"), 3))
print("rank by per-degree normalization (top 5):", r.sort_values("shared_down_per_100_outdegree_mean", ascending=False).TF.head(5).tolist())
for tf in ["STAT3", "CEBPB", "LYL1"]:
    x = g[g.TF == tf]; l, b = x.lung_delta.values, x.brain_delta.values
    print(tf, "down", int(((l < 0) & (b < 0)).sum()), "up", int(((l > 0) & (b > 0)).sum()),
          "down genes:", sorted(x.gene[(l < 0) & (b < 0)])[:40])

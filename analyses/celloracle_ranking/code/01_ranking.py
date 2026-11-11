# Advice section 4: full ranking of all common TFs from the frozen Entrez137 CellOracle outputs (no re-simulation).
# delta_X is in ln(1 + normalized count) space (normalize_per_cell -> CellOracle import log1p -> knn imputation),
# so delta = ln((1+x_sim)/(1+x_base)) is a signed log ratio; 100*delta/baseline is not an expression % change.
# Run from project root: _workspace/envs/celloracle/bin/python -I 03_candidate_analysis/CellOracle_ranking_scale_20261019/code/01_ranking.py
import pandas as pd, numpy as np
B = "03_candidate_analysis/CellOracle_AxisA_Entrez137_20261001/tables/"
O = "03_candidate_analysis/CellOracle_ranking_scale_20261019/"
p = pd.read_csv(B + "all_common_TF_paired_gene_responses.tsv.gz", sep="\t", low_memory=False)
p = p[p["mode"] == "main"]
p = p[p.in_model_lung.astype(str).eq("True") & p.in_model_brain.astype(str).eq("True")]
el = pd.read_csv(B + "all_TF_eligibility_and_execution.tsv.gz", sep="\t", low_memory=False)
deg = el.pivot_table(index="gene", columns="side", values="outgoing_nonzero", aggfunc="first")
rows = []
for tf, g in p.groupby("TF"):
    l, b = g.lung_delta.values, g.brain_delta.values
    rows.append(dict(TF=tf, n_evaluable=len(g),
        shared_down=int(((l < 0) & (b < 0)).sum()), shared_up=int(((l > 0) & (b > 0)).sum()),
        opposite=int(((l > 0) & (b < 0)).sum() + ((l < 0) & (b > 0)).sum()),
        lung_only_change=int(((l != 0) & (b == 0)).sum()), brain_only_change=int(((l == 0) & (b != 0)).sum()),
        mean_signed_delta_lung=l.mean(), mean_signed_delta_brain=b.mean(),
        mean_abs_delta_lung=np.abs(l).mean(), mean_abs_delta_brain=np.abs(b).mean(),
        degree_lung=deg.loc[tf, "lung"] if tf in deg.index else np.nan, degree_brain=deg.loc[tf, "brain"] if tf in deg.index else np.nan))
r = pd.DataFrame(rows)
r["shared_down_frac"] = r.shared_down / r.n_evaluable
r["shared_down_per_100_outdegree_mean"] = 100 * r.shared_down / r[["degree_lung", "degree_brain"]].mean(axis=1)
for s in ["lung", "brain"]:
    r[f"response_percentile_{s}"] = 100 * r[f"mean_abs_delta_{s}"].rank(pct=True)
    r[f"signed_percentile_{s}"] = 100 * r[f"mean_signed_delta_{s}"].rank(pct=True)
r = r.sort_values(["shared_down", "TF"], ascending=[False, True]).reset_index(drop=True)
r.insert(0, "rank_shared_down", r.shared_down.rank(ascending=False, method="min").astype(int))
assert len(r) == 357, len(r)
assert r.set_index("TF").loc["STAT3", "shared_down"] == 27 and r.set_index("TF").loc["FOS", "shared_down"] == 25
r.to_csv(O + "results/all357_TF_ranking.tsv", sep="\t", index=False)
# gene-level source data for redrawing 6d-f with signed log-ratio (natural log) and the invalid % for comparison
g = p[p.TF.isin(["STAT3", "CEBPB", "LYL1"])][["TF", "entrez_id", "gene", "lung_baseline", "lung_delta", "brain_baseline", "brain_delta"]].copy()
for s in ["lung", "brain"]:
    g[f"{s}_pct_of_log_baseline_old"] = 100 * g[f"{s}_delta"] / g[f"{s}_baseline"]
    g[f"{s}_linear_pct_change"] = 100 * (np.expm1(g[f"{s}_baseline"] + g[f"{s}_delta"]) - np.expm1(g[f"{s}_baseline"])) / np.expm1(g[f"{s}_baseline"])
g.to_csv(O + "results/Fig6def_source_signed_logratio.tsv", sep="\t", index=False)
cols = ["rank_shared_down", "TF", "n_evaluable", "shared_down", "shared_up", "opposite", "shared_down_frac", "degree_lung", "degree_brain", "shared_down_per_100_outdegree_mean", "response_percentile_lung", "response_percentile_brain"]
print(r[cols].head(12).round(2).to_string(index=False))
print("TFs with shared_down>0:", (r.shared_down > 0).sum(), "| shared_up>0:", (r.shared_up > 0).sum(), "| opposite>0:", (r.opposite > 0).sum())
print("Spearman shared_down vs mean outdegree:", round(r.shared_down.corr(r[["degree_lung", "degree_brain"]].mean(axis=1), method="spearman"), 3))
print("rank by per-degree normalization (top 10):", r.sort_values("shared_down_per_100_outdegree_mean", ascending=False).TF.head(10).tolist())
for tf in ["STAT3", "CEBPB", "LYL1"]:
    x = g[g.TF == tf]; nz = x[(x.lung_delta != 0) & (x.brain_delta != 0)]
    print(tf, "corr(old%, linear%) lung", round(nz.lung_pct_of_log_baseline_old.corr(nz.lung_linear_pct_change, method="spearman"), 3), "brain", round(nz.brain_pct_of_log_baseline_old.corr(nz.brain_linear_pct_change, method="spearman"), 3),
          "| median |old%|/|linear%| lung", round((nz.lung_pct_of_log_baseline_old.abs() / nz.lung_linear_pct_change.abs()).median(), 2))

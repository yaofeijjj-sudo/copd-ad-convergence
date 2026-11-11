from pathlib import Path
import hashlib, json, shutil
import numpy as np, pandas as pd

R = Path.cwd()
O = R / "results/gene_response_figures"
C = R / "results/virtual_knockout"
TFS = ["STAT3", "CEBPB", "LYL1"]
REUSED = []


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def build(tf, paired):
    a = paired[(paired.TF == tf) & (paired["mode"] == "main")].copy()
    a = a.dropna(subset=["lung_downstream_delta", "brain_downstream_delta"]).sort_values("gene")
    a = a[~a.forced_self_lung & ~a.forced_self_brain]
    assert not len(a[(a.lung_baseline <= 0) | (a.brain_baseline <= 0)]), "zero baseline"
    d = a[
        ["TF", "entrez_id", "gene", "gene_lung", "gene_brain", "ensembl_id_lung", "ensembl_id_brain"]
    ].copy()
    for side in ["lung", "brain"]:
        d[side + "_baseline"] = a[side + "_baseline"].values
        d[side + "_delta"] = a[side + "_downstream_delta"].values
        d[side + "_percent_of_baseline"] = 100 * d[side + "_delta"] / d[side + "_baseline"]
        d[side + "_rank"] = d[side + "_delta"].rank(method="average")
    d["pair_direction"] = np.select(
        [
            (d.lung_delta < 0) & (d.brain_delta < 0),
            (d.lung_delta > 0) & (d.brain_delta > 0),
            (d.lung_delta == 0) & (d.brain_delta == 0),
        ],
        ["shared_down", "shared_up", "both_zero"],
        default="other",
    )
    d["worse_organ_rank"] = d[["lung_rank", "brain_rank"]].max(axis=1)
    d["rank_sum"] = d.lung_rank + d.brain_rank
    d["label_order"] = 0
    sel = d[d.pair_direction == "shared_down"].sort_values(["worse_organ_rank", "rank_sum", "gene"]).head(6)
    d.loc[sel.index, "label_order"] = np.arange(1, len(sel) + 1)
    labels = [dict(TF=tf, gene=x.gene, kind="down") for _, x in sel.iterrows()]
    ups = d[d.pair_direction == "shared_up"].copy()
    ups["weak_percent"] = ups[["lung_percent_of_baseline", "brain_percent_of_baseline"]].min(axis=1)
    labels += [
        dict(TF=tf, gene=x.gene, kind="up")
        for _, x in ups.sort_values(["weak_percent", "gene"], ascending=[False, True]).head(2).iterrows()
    ]
    return d, labels


def lyl1_limits(d):
    lim = {"TF": "LYL1"}
    for side, lo, hi, ticks, neg, pos in [
        ("lung", "xmin", "xmax", "xticks", 5, 3),
        ("brain", "ymin", "ymax", "yticks", 1, 1),
    ]:
        v = d[side + "_percent_of_baseline"]
        scale = max(-v.min() / neg, v.max() / pos, 1e-9) * 1.03
        scale = np.ceil(scale / 3) * 3
        lim[lo], lim[hi] = -neg * scale, pos * scale
        step = 2 * scale if side == "lung" else scale / 2
        tt = np.arange(np.ceil(lim[lo] / step), np.floor(lim[hi] / step) + 1) * step
        lim[ticks] = ",".join(f"{x:g}" for x in tt)
    return lim


def run():
    paired_path = C / "tables/all_common_TF_paired_gene_responses.tsv.gz"
    paired = pd.read_csv(paired_path, sep="\t", dtype={"entrez_id": str})
    for name in ["source_data", "code", "qa"]:
        (O / name).mkdir(parents=True, exist_ok=True)
    limits_rel = pd.read_csv(R / "resources/gene_scatter_limits.tsv", sep="\t").set_index("TF")
    labels, limits, checks = [], [], {}
    for tf in TFS:
        d, lab = build(tf, paired)
        d.to_csv(O / f"source_data/{tf}_rank_rank_measured.tsv", sep="\t", index=False)
        lim = limits_rel.loc[tf].to_dict()
        lim["TF"] = tf
        checks[tf] = "donor means and relative changes recomputed"
        labels += lab
        limits.append(lim)
    pd.DataFrame(labels).to_csv(O / "source_data/scatter_selected_labels.tsv", sep="\t", index=False)
    cols = ["xmin", "xmax", "ymin", "ymax", "xticks", "yticks", "TF"]
    pd.DataFrame(limits)[cols].to_csv(O / "code/axis_limits.tsv", sep="\t", index=False)
    qa = dict(display_TFs=TFS, inputs={str(paired_path): sha(paired_path)}, checks=checks, limits=limits)
    (O / "qa/SOURCE_QA.json").write_text(json.dumps(qa, indent=2, ensure_ascii=False, default=str))
    print(json.dumps(qa, indent=1, ensure_ascii=False, default=str))


if __name__ == "__main__":
    run()

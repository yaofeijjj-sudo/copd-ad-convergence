from pathlib import Path
import pandas as pd, numpy as np, json

R = Path.cwd()
O = R / "results/virtual_knockout"


def run():
    (O / "source_data").mkdir(parents=True, exist_ok=True)
    response = pd.read_csv(
        O / "tables/all_AxisA137_gene_responses.tsv.gz", sep="\t", dtype={"entrez_id": str}
    )
    comp = pd.read_csv(O / "tables/lung_brain_full_comparison.tsv", sep="\t")
    common = sorted(comp.loc[comp.both_tested, "TF"])
    assert common
    identities = {}
    for side in ["lung", "brain"]:
        tab = pd.read_csv(O / f"tables/{side}_TF_eligibility.tsv", sep="\t")
        yes = tab[tab.KO_eligible]
        assert yes.annotation_symbol.is_unique
        identities[side] = yes.set_index("gene").annotation_symbol.to_dict()
    response["model_TF"] = response.TF
    response["TF"] = [identities[s][t] for s, t in zip(response.side, response.model_TF)]
    paired = []
    keys = ["TF", "mode", "entrez_id"]
    for mode in ["main", "randomized"]:
        frames = []
        for side in ["lung", "brain"]:
            a = (
                response[(response.side == side) & (response["mode"] == mode) & response.TF.isin(common)]
                .drop(columns="side")
                .copy()
            )
            assert not a.duplicated(keys).any()
            frames.append(a)
        z = frames[0].merge(
            frames[1], on=keys, how="outer", suffixes=("_lung", "_brain"), validate="one_to_one"
        )
        for name in ["current_symbol", "ensembl_id"]:
            if name == "current_symbol":
                assert z[name + "_lung"].equals(z[name + "_brain"])
        z["gene"] = z.current_symbol_lung
        for side in ["lung", "brain"]:
            for name in ["baseline", "delta", "downstream_delta", "simulated"]:
                z[side + "_" + name] = z[name + "_" + side]
        paired.append(z)
    paired = pd.concat(paired, ignore_index=True)
    assert len(paired) == len(common) * 137 * 2
    paired.to_csv(O / "tables/all_common_TF_paired_gene_responses.tsv.gz", sep="\t", index=False)
    rows = []
    for tf in common:
        a = paired[(paired.TF == tf) & (paired["mode"] == "main")]
        p = a.dropna(subset=["lung_downstream_delta", "brain_downstream_delta"])
        row = dict(
            TF=tf,
            common_measured_nonself=len(p),
            n_measured_genes=len(p),
            shared_down=int(((p.lung_downstream_delta < 0) & (p.brain_downstream_delta < 0)).sum()),
            shared_up=int(((p.lung_downstream_delta > 0) & (p.brain_downstream_delta > 0)).sum()),
        )
        for side in ["lung", "brain"]:
            b = a[a["in_model_" + side] & ~a["forced_self_" + side]].dropna(
                subset=[side + "_downstream_delta"]
            )
            ds = b[side + "_downstream_delta"]
            row.update(
                {
                    side + "_measured_nonself": len(b),
                    side + "_down": int((ds < 0).sum()),
                    side + "_mean_abs_delta": float(ds.abs().mean()),
                    side + "_n_nonzero": int((ds != 0).sum()),
                }
            )
        rows.append(row)
    breadth = pd.DataFrame(rows).sort_values("TF")
    nc = len(breadth)
    for side in ["lung", "brain"]:
        c = side + "_mean_abs_delta"
        rank = breadth[c].rank(method="min", ascending=False)
        rank = rank.mask(breadth[c] == 0, nc)
        breadth[side + "_strength_rank"] = rank
        breadth[side + "_strength_percentile"] = 100 * (nc - rank) / (nc - 1) if nc > 1 else 0
    breadth["shared_down_rank"] = breadth.shared_down.rank(method="min", ascending=False).astype(int)
    top = breadth[breadth.shared_down_rank <= 10].sort_values(["shared_down", "TF"], ascending=[False, True])
    breadth.to_csv(O / "source_data/all_common_TF_breadth_strength.tsv", sep="\t", index=False)
    top.to_csv(O / "source_data/overview.tsv", sep="\t", index=False)


if __name__ == "__main__":
    run()

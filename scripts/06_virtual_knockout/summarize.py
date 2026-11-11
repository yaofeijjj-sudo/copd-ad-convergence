from pathlib import Path
import pandas as pd

out = Path.cwd() / "results/virtual_knockout"
rows = []
identities = {}
for side in ["lung", "brain"]:
    axis = pd.read_csv(out / f"tables/{side}_AxisA_coverage.tsv", sep="\t", dtype={"entrez_id": str})
    ledger = pd.read_csv(out / f"tables/{side}_TF_eligibility.tsv", sep="\t")
    eligible = ledger.loc[ledger.KO_eligible]
    identities[side] = set(eligible.annotation_symbol)
    for tf in sorted(eligible.gene):
        for mode in ["main", "randomized"]:
            donor = pd.read_csv(out / f"tables/{side}_{tf}_{mode}_donor_genes.tsv.gz", sep="\t")
            means = donor.groupby("gene")[["baseline", "simulated", "delta"]].mean().reset_index()
            response = axis.merge(means, on="gene", how="left", validate="one_to_one")
            response["forced_self"] = response.gene.eq(tf)
            response["TF"] = tf
            response["side"] = side
            response["mode"] = mode
            response["downstream_delta"] = response.delta.mask(response.forced_self)
            assert response.loc[~response.in_model, "delta"].isna().all()
            rows.append(response)
pd.concat(rows, ignore_index=True).to_csv(
    out / "tables/all_AxisA137_gene_responses.tsv.gz", sep="\t", index=False
)
common = sorted(identities["lung"] & identities["brain"])
pd.DataFrame({"TF": common, "both_tested": True}).to_csv(
    out / "tables/lung_brain_full_comparison.tsv", sep="\t", index=False
)

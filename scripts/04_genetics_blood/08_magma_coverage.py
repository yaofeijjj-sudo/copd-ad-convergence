from pathlib import Path
import csv, json

R = Path.cwd()
O = R / "results/genetics_blood"
I = R / "results/module_overlap/results/cross"
rows = []
sets = {}
for disease in ["AD", "COPD"]:
    p = O / "results/magma/AD_total_N.genes.out" if disease == "AD" else O / "results/magma/COPD.genes.out"
    genes = {
        line.split()[0] for line in p.read_text().splitlines() if line.split() and line.split()[0].isdigit()
    }
    sets[disease] = {"_SET1_": genes, "_SET2_": genes}
for ax, tag in [("A", "_SET1_"), ("B", "_SET2_")]:
    for g in csv.DictReader((I / f"Axis{ax}_confirmed_genes.csv").open()):
        for disease in ["AD", "COPD"]:
            rows.append(
                dict(
                    disease=disease,
                    axis=f"Axis {ax}",
                    entrez_id=g["entrez_id"],
                    display_symbol=g["current_symbol"],
                    measured=g["entrez_id"] in sets[disease][tag],
                )
            )
with (O / "results/MAGMA_gene_coverage.csv").open("w") as f:
    w = csv.DictWriter(f, fieldnames=rows[0])
    w.writeheader()
    w.writerows(rows)
summary = []
for disease in ["AD", "COPD"]:
    for ax in ["Axis A", "Axis B"]:
        part = [r for r in rows if r["disease"] == disease and r["axis"] == ax]
        summary.append(
            dict(
                disease=disease,
                axis=ax,
                n_input=len(part),
                n_measured=sum(r["measured"] for r in part),
                missing_entrez=[r["entrez_id"] for r in part if not r["measured"]],
                missing_symbols=[r["display_symbol"] for r in part if not r["measured"]],
            )
        )
assert [x["n_measured"] for x in summary] == [130, 409, 130, 410], summary
(O / "qa/MAGMA_coverage_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))

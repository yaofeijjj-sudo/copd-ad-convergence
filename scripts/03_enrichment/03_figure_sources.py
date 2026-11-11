from pathlib import Path
import csv, collections

R = Path.cwd()
O = R / "results/enrichment"
D = O / "source_data"
D.mkdir(exist_ok=True)
a = list(csv.DictReader((O / "results/complete_family_all_terms.csv").open()))
lookup = {(r["gene_set"], r["source"], r["background_role"], r["term_id"]): r for r in a}
bg = R / "results/module_overlap/results/cross/shared_all_module_gene_ids.csv"
labels = {r["entrez_id"]: r["current_symbol"] for r in csv.DictReader(bg.open())}


def write(name, rows):
    with (D / name).open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


for axis in ["A", "B"]:
    terms = list(csv.DictReader((R / f"resources/enrichment_display_{axis}.csv").open()))
    members = {}
    updated = []
    for t in terms:
        db = {"GO BP": "GO_Biological_Process", "KEGG": "KEGGREST_hsa", "Hallmark": "MSigDB_Hallmark"}[
            t["database"]
        ]
        q = lookup[(f"Axis{axis}", db, "primary", t["term_id"])]
        qs = lookup[(f"Axis{axis}", db, "sensitivity", t["term_id"])]
        ids = q["overlap_genes"].split("/") if q["overlap_genes"] else []
        members[t["term_id"]] = ids
        t.update(
            gene_count=q["gene_count"],
            overlap_genes="/".join(labels[i] for i in ids),
            odds_ratio=q["odds_ratio"],
            fdr_bh=q["fdr_bh_complete_family"],
            neg_log10_fdr=q["neg_log10_fdr_complete_family"],
            gene_ratio=str(int(q["gene_count"]) / int(q["query_size"])),
            query_size=q["query_size"],
            universe_size=q["universe_size"],
            significant_primary=str(float(q["fdr_bh_complete_family"]) < 0.05).upper(),
            stability_class=(
                "significant_both"
                if float(q["fdr_bh_complete_family"]) < 0.05 and float(qs["fdr_bh_complete_family"]) < 0.05
                else "not_significant_both"
            ),
        )
        updated.append(t)
    degrees = collections.Counter(i for ids in members.values() for i in ids)
    ordered = sorted(degrees, key=lambda i: (-degrees[i], labels[i]))
    chosen = ordered[:18] if axis == "A" else [i for i in ordered if degrees[i] >= 2]
    chosen = sorted(chosen, key=lambda i: labels[i])
    genes = []
    links = []
    for n, i in enumerate(chosen, 1):
        ts = [t for t in terms if i in members[t["term_id"]]]
        rule = (
            "Top 18 by term recurrence; alphabetical tie-break"
            if axis == "A"
            else "All genes recurring in at least 2 terms"
        )
        genes.append(
            dict(
                axis=f"Axis {axis}",
                gene=labels[i],
                term_degree=degrees[i],
                mean_term_order=sum(int(t["display_order"]) for t in ts) / len(ts),
                selection_rule=rule,
                gene_order=n,
                entrez_id=i,
            )
        )
        for t in ts:
            links.append(
                dict(
                    display_order=t["display_order"],
                    axis=f"Axis {axis}",
                    database=t["database"],
                    term_id=t["term_id"],
                    term_name=t["term_name"],
                    gene=labels[i],
                    fdr_bh=t["fdr_bh"],
                    significant_primary=t["significant_primary"],
                    term_degree=degrees[i],
                    selection_rule=rule,
                    gene_order=n,
                    fdr_scope=t["fdr_scope"],
                    entrez_id=i,
                )
            )
    write(f"axis{axis}_functional_enrichment_complete_family_terms.csv", updated)
    write(f"axis{axis}_functional_enrichment_complete_family_links.csv", links)
    write(f"axis{axis}_functional_enrichment_selected_genes.csv", genes)

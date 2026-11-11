from pathlib import Path
import csv

root = Path.cwd()
out = root / "results/module_overlap_figures/source_data"
out.mkdir(parents=True, exist_ok=True)
with (root / "results/module_overlap/results/cross/all128_module_pairs_fisher_BH.csv").open() as f:
    rows = list(csv.DictReader(f))
for name, selected in [
    ("panel_a_all128.csv", rows),
    (
        "panels_bc_venn_counts.csv",
        [
            x
            for x in rows
            if (x["lung_module"], x["brain_module"]) in [("red", "blue"), ("turquoise", "turquoise")]
        ],
    ),
]:
    with (out / name).open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(selected[0]))
        writer.writeheader()
        writer.writerows(selected)

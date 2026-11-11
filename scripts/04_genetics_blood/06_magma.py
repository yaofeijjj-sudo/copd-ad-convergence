from pathlib import Path
import csv
import os
import subprocess

root = Path.cwd()
out = root / "results/genetics_blood"
inputs = out / "results/magma_inputs"
results = out / "results/magma"
for p in [inputs, results, out / "logs"]:
    p.mkdir(parents=True, exist_ok=True)
binary = root / "data/tools/magma/magma"
ld = root / "data/tools/magma/g1000_eur/g1000_eur"
genes = root / "data/tools/magma/NCBI37/NCBI37.3.gene.loc"
annotation = inputs / "annotation.genic_only"
env = os.environ.copy()
env["LC_ALL"] = "C"


def run(arguments, name):
    with (out / "logs" / f"{name}.log").open("w") as log:
        subprocess.run(
            [str(binary), *map(str, arguments)], stdout=log, stderr=subprocess.STDOUT, check=True, env=env
        )


if not Path(str(annotation) + ".genes.annot").exists():
    run(["--annotate", "--snp-loc", str(ld) + ".bim", "--gene-loc", genes, "--out", annotation], "annotation")
sets = inputs / "axes.sets"
with sets.open("w") as f:
    for axis in ["A", "B"]:
        with (root / f"results/module_overlap/results/cross/Axis{axis}_confirmed_genes.csv").open() as src:
            ids = [x["entrez_id"] for x in csv.DictReader(src)]
        f.write("Axis" + axis + "\t" + "\t".join(ids) + "\n")
for disease in ["AD", "COPD"]:
    prefix = results / ("AD_total_N" if disease == "AD" else "COPD")
    if not Path(str(prefix) + ".genes.raw").exists():
        run(
            [
                "--bfile",
                ld,
                "--gene-annot",
                str(annotation) + ".genes.annot",
                "--pval",
                inputs / f"{disease}.magma.pval",
                "ncol=N",
                "--out",
                prefix,
            ],
            disease + "_gene",
        )
    run(
        [
            "--gene-results",
            str(prefix) + ".genes.raw",
            "--set-annot",
            sets,
            "--out",
            results / f"{disease}_axes",
        ],
        disease + "_gsa",
    )

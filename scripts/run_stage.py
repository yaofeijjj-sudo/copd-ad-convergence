from pathlib import Path
import argparse, json, os, subprocess, sys


repo = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument("stage")
p.add_argument("--config", type=Path, default=repo / "config/paths.yml")
a = p.parse_args()
paths = dict(
    line.split(":", 1)
    for line in a.config.read_text().splitlines()
    if line.strip() and not line.startswith("#")
)
paths = {k: v.strip().strip(chr(34)).strip(chr(39)) for k, v in paths.items()}


def resolve(name):
    x = Path(paths[name]).expanduser()
    return x if x.is_absolute() else repo / x


work = resolve("output_root")
work.mkdir(parents=True, exist_ok=True)
data = resolve("data_root")
data.mkdir(parents=True, exist_ok=True)
for name, target in [
    ("data", data),
    ("scripts", repo / "scripts"),
    ("R", repo / "R"),
    ("config", repo / "config"),
    ("resources", repo / "resources"),
    ("validation", repo / "validation"),
]:
    link = work / name
    if not link.exists():
        link.symlink_to(target, target_is_directory=True)
for name in [
    "tissue_scope/metadata",
    "wgcna/lung",
    "wgcna/brain",
    "module_overlap/results/identity",
    "module_overlap/results/cross",
    "module_overlap/provenance",
    "module_overlap/logs",
    "module_overlap/qa",
    "enrichment/results",
    "enrichment/provenance",
    "genetics_blood/results/magma_inputs",
    "genetics_blood/results/magma",
    "genetics_blood/provenance",
    "genetics_blood/logs",
    "genetics_blood/qa",
    "single_cell/results",
    "single_cell/provenance",
    "single_cell/qa",
    "virtual_knockout/objects",
    "virtual_knockout/tables",
    "virtual_knockout/protocol",
    "virtual_knockout/sources",
    "virtual_knockout/qa",
    "virtual_knockout/logs",
    "cache",
]:
    (work / "results" / name).mkdir(parents=True, exist_ok=True)
stages = json.loads((repo / "config/stages.json").read_text())
env = os.environ.copy()
env["WGCNA_AD_SKIP_RESIDUALIZATION"] = "1"
env["WGCNA_SKIP_PRESERVATION"] = "1"
for step in stages[a.stage]:
    program = "Rscript" if step[0].endswith(".R") else sys.executable
    subprocess.run([program, *step], cwd=work, env=env, check=True)

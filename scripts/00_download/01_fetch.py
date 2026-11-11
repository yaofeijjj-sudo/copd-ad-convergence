"""Download the public inputs listed in data/manifest.tsv into data_root and check them against the recorded sha256.

Entries whose source is a URL are downloaded. GSE302339 count matrices are extracted from the RAW tar, the NCBI Gene
snapshot is copied from resources/, and the Ensembl 116 transcript table is derived from the GTF. Platform tables,
KEGG and Hallmark annotations are written by 02_annotations.R. MAGMA and its reference files are listed at the
end for manual download.
"""

import csv, gzip, hashlib, shutil, sys, tarfile, urllib.request
from pathlib import Path

repo = Path(__file__).resolve().parents[2]
paths = dict(
    line.split(":", 1) for line in (repo / "config/paths.yml").read_text().splitlines() if ":" in line
)
root = Path(paths["data_root"].strip().strip("\"'"))
root = root if root.is_absolute() else repo / root
rows = list(csv.DictReader(open(repo / "data/manifest.tsv"), delimiter="\t"))
GTF = "https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def download(url, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(dst.suffix + ".part")
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f, 1 << 24)
    tmp.rename(dst)


def tx2gene(dst):
    tmp = root / "reference/Homo_sapiens.GRCh38.116.gtf.gz"
    if not tmp.exists():
        download(GTF, tmp)
    with gzip.open(tmp, "rt") as f, gzip.open(dst, "wt") as out:
        for line in f:
            x = line.split("\t")
            if len(x) > 8 and x[2] == "transcript":
                a = dict(kv.strip().split(" ", 1) for kv in x[8].split(";") if kv.strip())
                out.write(a["transcript_id"].strip('"') + "\t" + a["gene_id"].strip('"') + "\n")


manual = []
for r in rows:
    dst = root / r["path"]
    src = r["source"]
    if dst.exists():
        continue
    if src.startswith(("http://", "https://", "ftp://")):
        print("download", r["path"])
        download(src.replace("ftp://", "https://", 1), dst)
    elif src.startswith("resources/"):
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / src.split(" ")[0], dst)
    elif src == "extracted from GSE302339_RAW.tar":
        with tarfile.open(root / "raw/single_cell/copd_lung_gse302339/GSE302339_RAW.tar") as t:
            t.extract(Path(r["path"]).name, dst.parent)
    elif r["path"].endswith("ensembl116_tx2gene.tsv.gz"):
        tx2gene(dst)
    elif not src.startswith(("GEOquery", "KEGGREST", "msigdbr", "submitter platform table")):
        manual.append(r)

# Derived gzip files are compared on their content, since gzip headers carry timestamps.
check = []
for r in rows:
    p = root / r["path"]
    if not p.exists():
        check.append((r["path"], "missing"))
        continue
    if sha256(p) == r["sha256"]:
        check.append((r["path"], "ok"))
    elif r["path"].endswith(".gz") and not r["source"].startswith(("http", "ftp")):
        check.append((r["path"], "content differs from recorded gzip; compare decompressed"))
    else:
        check.append((r["path"], "sha256 differs"))
with open(root / "data_check.tsv", "w") as f:
    f.write("path\tstatus\n")
    f.writelines(f"{a}\t{b}\n" for a, b in check)
bad = [c for c in check if c[1] != "ok"]
print(len(rows) - len(bad), "of", len(rows), "files match the manifest; details in", root / "data_check.tsv")
for r in manual:
    print("manual:", r["path"], "<-", r["source"])
sys.exit(1 if bad else 0)

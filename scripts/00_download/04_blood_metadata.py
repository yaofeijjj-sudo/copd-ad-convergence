"""Resolve GEO sample IDs against the author's published clinical table."""

import csv
import gzip
import hashlib
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

ROOT = Path.cwd()
RAW = ROOT / "data/raw/bulk_blood_external/GSE248417"
OUT = ROOT / "results/blood_metadata"
OUT.mkdir(parents=True, exist_ok=True)
MATRIX = RAW / "GSE248417_ADvsControl.mRNA.gene_level.deg.txt.gz"
SERIES = RAW / "GSE248417_series_matrix.txt.gz"
CLINICAL = RAW / "supplementary_clinical_data.docx"
AXIS = ROOT / "resources/AxisA_confirmed_genes.csv"


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


with gzip.open(SERIES, "rt", encoding="utf-8") as f:
    geo_rows = list(csv.reader(f, delimiter="\t"))
geo = {}
for row in geo_rows:
    if row and row[0].startswith("!Sample_"):
        geo.setdefault(row[0], []).append([v.strip('"') for v in row[1:]])

with ZipFile(CLINICAL) as z:
    root = ET.fromstring(z.read("word/document.xml"))
ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
table = root.find(".//w:tbl", ns)
assert table is not None
clinical_rows = []
for tr in table.findall("./w:tr", ns):
    clinical_rows.append(
        [
            " ".join((x.text or "") for x in tc.findall(".//w:t", ns)).strip()
            for tc in tr.findall("./w:tc", ns)
        ]
    )
headers = clinical_rows[0]
clinical = {r["Sample"]: r for r in (dict(zip(headers, row)) for row in clinical_rows[1:])}
assert len(clinical) == 100

descriptions = geo["!Sample_description"][0]
gsm = geo["!Sample_geo_accession"][0]
disease = [x.split(":", 1)[1].strip() for x in geo["!Sample_characteristics_ch1"][1]]
assert len(descriptions) == len(gsm) == len(disease) == 98
assert len(set(descriptions)) == 98

with gzip.open(MATRIX, "rt", encoding="utf-8") as f:
    matrix_header = next(csv.reader(f, delimiter="\t"))
sample_cols = matrix_header[4:-5]
assert len(sample_cols) == 98 and set(sample_cols) == set(descriptions)

records = []
for sample_id, accession, label in zip(descriptions, gsm, disease):
    assert label in {"AD", "Control"}
    source_id = sample_id.replace("AD", "A", 1) if sample_id.startswith("AD") else sample_id
    metadata = clinical.get(source_id)
    if metadata is not None:
        assert metadata["Sample Type"] == label
        assert metadata["Age"].isdigit() and metadata["Sex"] in {"F", "M"}
    include_primary = metadata is not None and "_" not in sample_id
    records.append(
        {
            "geo_accession": accession,
            "sample_id": sample_id,
            "supplement_id": source_id if metadata else "",
            "group": label,
            "age": metadata["Age"] if metadata else "",
            "sex": metadata["Sex"] if metadata else "",
            "include_primary": str(include_primary).lower(),
            "exclusion_reason": "" if include_primary else "composite_unmatched_control_id",
        }
    )

assert sum(r["group"] == "AD" for r in records) == 49
assert sum(r["group"] == "Control" for r in records) == 49
assert sum(r["include_primary"] == "true" for r in records) == 97
assert [r["sample_id"] for r in records if r["include_primary"] == "false"] == ["C8202_7603"]

with (OUT / "sample_manifest_gse248417.csv").open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=records[0].keys())
    writer.writeheader()
    writer.writerows(records)

with (OUT / "INPUT_SHA256.tsv").open("w", encoding="utf-8") as f:
    f.write("path\tsha256\n")
    for path in (MATRIX, SERIES, CLINICAL, AXIS):
        f.write(f"{path.relative_to(ROOT)}\t{sha256(path)}\n")

print(f"GEO: {len(records)} samples; primary: {sum(r['include_primary'] == 'true' for r in records)}")
print("Unmatched clinical ID:", [r["sample_id"] for r in records if r["include_primary"] == "false"])
print("Clinical-only IDs:", sorted(set(clinical) - {r["supplement_id"] for r in records}))

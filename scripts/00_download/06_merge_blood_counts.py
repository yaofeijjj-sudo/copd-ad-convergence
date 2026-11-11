from pathlib import Path
import pandas as pd

root = Path.cwd()
raw = root / "data/raw/bulk_blood_external/GSE248493"
out = root / "results/blood_ranks/inputs"
out.mkdir(parents=True, exist_ok=True)
parts = []
manifest = []
for file in sorted(raw.glob("GSE248493_count_*.xlsx")):
    d = pd.read_excel(file, index_col=0)
    counts = d[[c for c in d if c.startswith("DK")]].copy()
    assert counts.index.is_unique and not counts.isna().any().any()
    if parts:
        assert parts[0].index.equals(counts.index)
    parts.append(counts)
    manifest.extend({"source_file": file.name, "count_sample_id": c} for c in counts)
assert len(parts) == 4
counts = pd.concat(parts, axis=1)
assert counts.columns.is_unique
counts.index.name = "gene_symbol"
counts.to_csv(out / "GSE248493_merged_counts.csv.gz")
pd.DataFrame(manifest).to_csv(out / "GSE248493_count_file_manifest.csv", index=False)

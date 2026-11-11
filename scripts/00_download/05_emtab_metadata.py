#!/usr/bin/env python3
"""Freeze subject inclusion from ArrayExpress SDRF and processed-matrix header."""
import csv
from pathlib import Path
from collections import Counter

ROOT = Path.cwd()
RAW = ROOT / "data/raw/bulk_blood_external/E-MTAB-5278"
OUT = ROOT / "results/blood_metadata"
OUT.mkdir(parents=True, exist_ok=True)
SDRF = RAW / "E-MTAB-5278.sdrf.txt"
MATRIX = RAW / "H1_QASMC_Blood_Hs_MRNA_normMatrix.txt"

with MATRIX.open() as stream:
    matrix_columns = stream.readline().rstrip("\n\r").split("\t")[1:]
assert len(matrix_columns) == len(set(matrix_columns)) == 224
with SDRF.open() as stream:
    source = list(csv.DictReader(stream, delimiter="\t"))
assert len(source) == 234
assert Counter(row["Characteristics[disease]"] for row in source) == {
    "normal": 179,
    "chronic obstructive pulmonary disease": 55,
}

rows = []
for row in source:
    disease = row["Characteristics[disease]"]
    smoking = row["Characteristics[clinical history]"]
    has_expression = row["Array Data File"] in matrix_columns
    is_current_smoker = smoking == "smoker"
    include = has_expression and is_current_smoker
    rows.append(
        {
            "donor_id": row["Source Name"],
            "external_sample_id": row["Comment[external sample id]"],
            "array_file": row["Array Data File"],
            "group": "Case" if disease == "chronic obstructive pulmonary disease" else "Control",
            "age": row["Characteristics[age]"],
            "sex": "M" if row["Characteristics[sex]"] == "male" else "F",
            "smoking_status": smoking,
            "race": row["Characteristics[ethnic group]"],
            "present_in_processed_matrix": str(has_expression).lower(),
            "include_primary": str(include).lower(),
            "exclusion_reason": (
                "included"
                if include
                else (
                    "not_current_smoker"
                    if not is_current_smoker
                    else "absent_from_submitted_processed_matrix"
                )
            ),
        }
    )

selected = [row for row in rows if row["include_primary"] == "true"]
assert Counter(row["group"] for row in selected) == {"Case": 51, "Control": 58}
assert len({row["donor_id"] for row in selected}) == 109
assert len({row["array_file"] for row in selected}) == 109
assert all(row["age"].isdigit() and row["sex"] in {"M", "F"} for row in selected)
assert len({row["sex"] for row in selected}) == 2
with (OUT / "sample_manifest_emtab5278.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
print("Selected 51 COPD and 58 smoker controls; 109 unique donors")

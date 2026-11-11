"""2026-10-29: write the SEA-AD region as DLPFC in S6 Table, matching the manuscript.

The frozen SEA-AD result files label the region with the Allen ontology code "DFC" (cohort labels such as
"Brain (SEA-AD DFC)" and the CellOracle input description). Only displayed text is changed; values and the
CPS_Local variable name are untouched. Runs after append_S6g_alveolar_macrophage.py (see run_all.sh). Idempotent.
Run from the project root:
    _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/update_S6_region_label.py
"""
import re
from pathlib import Path

import openpyxl

P = Path.cwd() / "06_resubmission_20261008/code/supplement/parts/S6_Table.xlsx"
wb = openpyxl.load_workbook(P)
n = 0
for ws in wb:
    for row in ws.iter_rows():
        for c in row:
            if isinstance(c.value, str) and re.search(r"\bDFC\b", c.value):
                c.value = re.sub(r"\bDFC\b", "DLPFC", c.value)
                n += 1
wb.save(P)
print(f"S6 Table: {n} cells relabelled DFC -> DLPFC")

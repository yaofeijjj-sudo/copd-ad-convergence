"""Update main-figure citations inside the S-table parts after the 2026-11-05 renumbering (Fig 1-6 became Fig 2-8:
new Fig 1 = study overview; former Fig 3 split into Fig 4 and Fig 5). Runs after update_parts_20261105.py. Idempotent.
Run from the project root:
    _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/update_figure_numbers_20261105.py
"""
from pathlib import Path
import openpyxl

P = Path.cwd() / "06_resubmission_20261008/code/supplement/parts"
REPL = [("(as in Fig 1a and 1b)", "(as in Fig 2a and 2b)"), ("(metafor 5.0.1; Fig 2b)", "(metafor 5.0.1; Fig 3c)"),
        ("Display names follow Fig 3;", "Display names follow Fig 4;"), ("plus the Fig 1e covariates", "plus the Fig 2e covariates"),
        ("covariates of Fig 1e", "covariates of Fig 2e"), ("a–e; Fig 4c)", "a–e; Fig 6c)"), ("(UCell, Fig localization)", "(UCell, Fig 4)")]
done = {a: 0 for a, _ in REPL}
for f in sorted(P.glob("S*_Table.xlsx")):
    wb = openpyxl.load_workbook(f); changed = False
    for ws in wb.worksheets:
        for row in ws.iter_rows():
            for c in row:
                if isinstance(c.value, str):
                    v = c.value
                    for a, b in REPL:
                        if a in v:
                            v = v.replace(a, b); done[a] += 1
                    if v != c.value:
                        c.value = v; changed = True
    if changed:
        wb.save(f)
print(done)

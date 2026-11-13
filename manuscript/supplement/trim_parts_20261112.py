"""11-12: legends state what was done (no 'not analyzed' wording). Keep in S2 File only what the manuscript, its figures or the S1 File methods report.

Author's rule (2026-11-12): Methods and supplements only need to let readers reproduce the reported results; checks and
sensitivity analyses that no reported number depends on are removed. This script edits the table parts in place, after
the other update scripts of run_all.sh and before build_S1_S2_File.py. It is idempotent.

Removed sheets: S1b, S1c, S3b, S4b, S4c, S4d, S5b (member effects are in S1 Data, Fig3d), S5c, S6c, S6d, S6e, S7c, S7e,
S7g, S10f, S11b. Sub-parts are re-lettered (S4e->S4b, S6f->S6c, S6g->S6d, S7d->S7c, S7f->S7d).
Removed rows/columns: chronic hypersensitivity pneumonitis, confounded with batch and not analyzed (S4a); modified Hartung-Knapp (S5a); FEV1 and three-group models (S6b); alveolar macrophage
sensitivity rows and Hodges-Lehmann columns (S6g); reference version A and 'all proportions except reference' (S7a, S7b,
S7d); dtangle and multi-type rows (S7f); shared-genes-removed rows and VIF (S9b, S9c); set-only APOE variant (S10a);
all-1,000 BRETIGEA set and the descriptive non-reference rows (S10b, S10c, S10d); expression-only P (S9d); the out-degree correlation sentence (S12a; out-degree itself is kept because Fig 8b shows it); non-primary models, the prior
without Axis A targets and the all-TF FDR (S13a).
2026-11-13 (author): the Results report the fgsea test of the blood datasets again, so the fgsea columns of S11a are
kept and a BH FDR across the six datasets is added next to them, as for CAMERA.

Run from the project root:
  _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/trim_parts_20261112.py
"""
import copy
from pathlib import Path

import openpyxl

PARTS = Path(__file__).resolve().parent / "parts"
HEADER_ROW, FIRST = 4, 5


def header(ws):
    return [ws.cell(HEADER_ROW, c).value for c in range(1, ws.max_column + 1)]


def keep_rows(ws, col_name, keep):
    c = header(ws).index(col_name) + 1
    for r in range(ws.max_row, FIRST - 1, -1):
        v = ws.cell(r, c).value
        if v is not None and not keep(str(v)):
            ws.delete_rows(r)


def drop_cols(ws, pred):
    keep = {a: (ws[a].value, copy.copy(ws[a].font), copy.copy(ws[a].alignment)) for a in ("A1", "A2", "A3")}
    merged = [m for m in ws.merged_cells.ranges if m.min_row == 2]
    for m in merged:
        ws.unmerge_cells(str(m))
    h = header(ws)
    for c in range(len(h), 0, -1):
        if h[c - 1] is not None and pred(c - 1, h):
            ws.delete_cols(c)
    for a, (v, f, al) in keep.items():  # title, legend and spacer live in column A; keep them if column A was deleted
        ws[a].value, ws[a].font, ws[a].alignment = v, f, al
    if merged:
        ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=ws.max_column)


def text(ws, cell, old, new):
    v = ws[cell].value
    if old in v:
        ws[cell].value = v.replace(old, new)


def edit(n, fn):
    p = PARTS / f"S{n}_Table.xlsx"
    wb = openpyxl.load_workbook(p)
    if "trimmed-20261112" in (wb.properties.keywords or ""):  # already trimmed: running again changes nothing
        return
    fn(wb)
    wb.properties.keywords = ((wb.properties.keywords or "") + " trimmed-20261112").strip()
    wb.save(p)


def remove(wb, *names):
    for s in names:
        if s in wb.sheetnames:
            del wb[s]


def relabel(wb, old, new, letter):
    if old in wb.sheetnames:
        ws = wb[old]
        ws.title = new
        ws["A1"].value = ws["A1"].value.replace(f"Table ({old[-1]})", f"Table ({letter})")


def s1(wb):
    remove(wb, "S1b", "S1c")
    text(wb["S1a"], "A1", ", and sensitivity analyses with the grey-module background and with restored age outliers.", ".")


def s3(wb):
    remove(wb, "S3b")
    ws = wb["S3a"]
    text(ws, "A2", "(a) Set level. ", "")
    text(ws, "A2", "computed without disease labels", "computed across all samples")
    text(ws, "A2", ", which exclude GSM1176257 from GSE48350 (S4 Table (b))", "")
    text(ws, "A2", " (b) Member level: each member's median Spearman correlation with the other members in each dataset. "
         "'consistent' = median > 0 in every dataset where measured and measured in >= 6 of 8 datasets; labels are "
         "descriptive and do not redefine Axis A.", "")


def s4(wb):
    remove(wb, "S4b", "S4c", "S4d")
    relabel(wb, "S4e", "S4b", "b")
    text(wb["S4b"], "A2", "P values are unadjusted, and these comparisons are outside the BH-corrected disease scope families.",
         "P values are unadjusted.")
    ws = wb["S4a"]
    text(ws, "A1", ", data processing and donor checks for", " and data processing for")
    text(ws, "A2", " and its sex-gene expression is male, consistent with that control (S4 Table b)", "")
    keep_rows(ws, "Contrast", lambda v: "hypersensitivity" not in v.lower())
    text(ws, "A2", " Chronic hypersensitivity pneumonitis is shown descriptively without FDR: in GSE150910 all its samples were "
         "processed on dates without controls, so disease and processing batch cannot be separated.", "")


def s5(wb):
    remove(wb, "S5b", "S5c")
    ws = wb["S5a"]
    text(ws, "A1", "Sensitivity analyses of Axis A replication and member-level disease effects.",
         "Sensitivity analyses of the pooled Axis A effect in AD.")
    text(ws, "A2", " REML + modified Hartung-Knapp (mKH): the Hartung-Knapp variance scaling factor truncated at 1 "
         "(Röver et al. 2015); here the factor was 1.006, so the interval equals the unmodified Hartung-Knapp interval.", "")
    keep_rows(ws, "Analysis", lambda v: "modified" not in v)


def s6(wb):
    remove(wb, "S6c", "S6d", "S6e")
    ws = wb["S6a"]
    text(ws, "A1", "(UCell, pseudobulk, alveolar macrophages, emphysema subtypes and lung function) and quality control and "
         "annotation of the lung single-cell data.",
         "(UCell, pseudobulk and alveolar macrophages) and annotation of the lung single-cell data.")
    text(ws, "A2", " (d) Per-library quality control of the lung data. (e) Lung myeloid test after removing flagged "
         "libraries. (f) CellTypist labels and broad classes.", " (c) CellTypist labels and broad classes.")
    ws = wb["S6b"]
    keep_rows(ws, "Model", lambda v: "FEV1" not in v and "3-group" not in v)
    text(ws, "A2", " E-COPD/NE-COPD: emphysematous and non-emphysematous COPD as defined by the original authors;", "")
    ws = wb["S6g"]
    keep_rows(ws, "Analysis", lambda v: v == "Alveolar macrophages (primary)")
    hl = [i for i, h in enumerate(header(ws)) if h and h.startswith("Hodges")]
    if hl:
        drop_cols(ws, lambda i, h: hl[0] <= i <= hl[0] + 2)
    text(ws, "A2", " with the Hodges–Lehmann location shift and its 95% CI", "")
    text(ws, "A2", " (50 in the stated sensitivity analysis). Flagged libraries: S6 Table (d). Alveolar macrophage family: "
         "'Alveolar macrophages', 'Alveolar Mph CCL3+', 'Alveolar Mph proliferating' and 'Alveolar Mph MT-positive'. "
         "Alveolar macrophages as a share of myeloid cells per donor, median (IQR): COPD 0.49 (0.39–0.61), control 0.48 "
         "(0.29–0.66).", ".")
    text(wb["S6f"], "A2", "best-match mode, no majority voting)", "best-match mode)")
    text(wb["S6f"], "A2", " Other/Progenitor (hematopoietic stem cells) was not analyzed as a class.", "")
    relabel(wb, "S6f", "S6c", "c")
    relabel(wb, "S6g", "S6d", "d")


def s7(wb):
    remove(wb, "S7c", "S7e", "S7g")
    for s in ["S7a", "S7b", "S7d"]:
        ws = wb[s]
        keep_rows(ws, "Reference version", lambda v: v == "B")
        if "Model" in header(ws):
            keep_rows(ws, "Model", lambda v: not v.startswith("all proportions"))
        drop_cols(ws, lambda i, h: h[i] == "Reference version")
    text(wb["S7a"], "A2", "Version B (primary): signature built after removing all 137 Axis A genes; version A: all genes. "
         "'target proportion' adds the standardized myeloid (lung) or immune (brain) proportion to the replication model; "
         "'all proportions except reference' adds every estimated proportion except the named reference class.",
         "The CIBERSORTx signature matrices were built after removing all 137 Axis A genes. 'target proportion' adds the "
         "standardized myeloid (lung) or immune (brain) proportion to the replication model.")
    ws = wb["S7f"]
    keep_rows(ws, "Role", lambda v: v == "primary")
    keep_rows(ws, "Method", lambda v: v == "BRETIGEA")
    drop_cols(ws, lambda i, h: h[i] in ("Role", "Marker version"))
    text(ws, "A1", "Discovery cohorts: module–disease effects adjusted for proportions estimated with dtangle (lung) and "
         "BRETIGEA (brain).", "Brain discovery cohort GSE33000: blue module–AD effect adjusted for the BRETIGEA microglial "
         "proportion.")
    relabel(wb, "S7d", "S7c", "c")
    relabel(wb, "S7f", "S7d", "d")


def s9(wb):
    for s in ["S9b", "S9c"]:
        ws = wb[s]
        keep_rows(ws, "Analysis", lambda v: "shared genes removed" not in v)
        drop_cols(ws, lambda i, h: h[i] in ("VIF", "Maximum VIF"))
    text(wb["S9b"], "A2", "; 'shared genes removed from both': the same after removing shared members from both scores. "
         "VIF, variance inflation factor of the signature term;", ";")
    text(wb["S9c"], "A2", "; maximum VIF", "")
    drop_cols(wb["S9d"], lambda i, h: h[i] == "Empirical P, matched on expression only")
    text(wb["S9d"], "A2", "not smaller than", "greater than or equal to")


def s10(wb):
    remove(wb, "S10f")
    ws = wb["S10a"]
    keep_rows(ws, "Variant", lambda v: v != "V1_set_only")
    text(ws, "A1", ", restricted background and annotation windows)", " and restricted background)")
    text(ws, "A2", "'set_only': region genes removed from the gene sets, background unchanged; 'whole_region': the 69 region "
         "genes removed from the whole analysis with MAGMA --settings gene-exclude (reported in the main text).",
         "'whole_region': the 69 region genes removed from the whole analysis with MAGMA --settings gene-exclude.")
    text(ws, "A2", "the top 250 (primary) or all 1,000 (sensitivity) human microglia markers", "the top 250 human microglia markers")
    for s in ["S10b", "S10c"]:
        keep_rows(wb[s], "Reference set", lambda v: v == "top250")
    keep_rows(wb["S10b"], "Variable", lambda v: v != "AxisA_nonRef")
    keep_rows(wb["S10d"], "Variable", lambda v: v != "AxisA_nonPatir")
    text(wb["S10b"], "A2", " AxisA_nonRef = Axis A members outside the reference set (descriptive).", "")


def s11(wb):
    remove(wb, "S11b")
    ws = wb["S11a"]
    text(ws, "A1", ", participant-level effects and duplicate-donor checks for", " and participant-level effects for")
    text(ws, "A1", "CAMERA tests and participant-level", "fgsea and CAMERA tests and participant-level")
    text(ws, "A2", " fgsea (fgseaMultilevel on moderated t) NES is shown for direction only; fgsea treats genes as "
         "independent, so its P values are smaller than those of CAMERA, which accounts for inter-gene correlation.",
         " fgsea: fgseaMultilevel on genes ranked by moderated t statistics; BH across the six datasets.")
    h = header(ws)
    rows_ = [r for r in range(FIRST, ws.max_row + 1) if ws.cell(r, 1).value is not None]

    def bh(col):
        p = [float(ws.cell(r, h.index(col) + 1).value) for r in rows_]
        out, m = [None] * len(p), 1.0
        for rank, i in sorted(enumerate(sorted(range(len(p)), key=lambda k: p[k]), 1), reverse=True):
            m = min(m, p[i] * len(p) / rank)
            out[i] = m
        return out

    assert len(rows_) == 6 and h[-1] == "fgsea P", (len(rows_), h)
    assert all(abs(a - float(ws.cell(r, h.index("CAMERA BH FDR (6 datasets)") + 1).value)) < 1e-12
               for a, r in zip(bh("CAMERA P"), rows_))  # same BH as the CAMERA column
    merged = [m for m in ws.merged_cells.ranges if m.min_row == 2]
    for m in merged:
        ws.unmerge_cells(str(m))
    src, new = len(h), len(h) + 1
    for r, v in zip([HEADER_ROW] + rows_, ["fgsea BH FDR (6 datasets)"] + bh("fgsea P")):
        a, b = ws.cell(r, src), ws.cell(r, new)
        b.value, b.font, b.alignment, b.border, b.fill, b.number_format = (
            v, copy.copy(a.font), copy.copy(a.alignment), copy.copy(a.border), copy.copy(a.fill), a.number_format)
    ws.column_dimensions[openpyxl.utils.get_column_letter(new)].width = ws.column_dimensions[openpyxl.utils.get_column_letter(src)].width
    if merged:
        ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=ws.max_column)


def s12(wb):
    ws = wb["S12a"]
    drop_cols(ws, lambda i, h: h[i] == "Shared decreases per 100 out-degree")
    text(ws, "A2", " Across the 357 TFs, the number of shared decreases increases with the total out-degree (lung plus "
         "brain; Spearman ρ = 0.68).", "")


def s13(wb):
    ws = wb["S13a"]
    keep_rows(ws, "Model", lambda v: v.startswith("primary"))
    prior = [h for h in header(ws) if h == "Regulatory prior"][0]
    full = ws.cell(FIRST, header(ws).index(prior) + 1).value
    vals = {str(ws.cell(r, header(ws).index(prior) + 1).value) for r in range(FIRST, ws.max_row + 1)}
    keep_rows(ws, prior, lambda v: "remov" not in v.lower() and "without" not in v.lower())
    drop_cols(ws, lambda i, h: h[i] == "BH FDR, all TFs")
    text(ws, "A2", ", with the full CollecTRI prior or the prior after removing Axis A targets", " with the CollecTRI prior")
    text(ws, "A2", "design matrix of the pseudobulk tests (S6 Table b)", "design matrix of the primary pseudobulk tests (S6 Table b)")
    text(ws, "A2", " (primary family) and within all TFs", "")
    print("S13 prior values seen:", sorted(vals), "first:", full)


for n, fn in [(1, s1), (3, s3), (4, s4), (5, s5), (6, s6), (7, s7), (9, s9), (10, s10), (11, s11), (12, s12), (13, s13)]:
    edit(n, fn)
for n in range(1, 14):
    wb = openpyxl.load_workbook(PARTS / f"S{n}_Table.xlsx")
    print(f"S{n}", [(ws.title, ws.max_row - HEADER_ROW, ws.max_column) for ws in wb.worksheets])

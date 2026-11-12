#!/bin/sh
# Rebuild all supporting information (S1 File, S2 File, S1 Data) from the frozen result files, in dependency order.
# Run from the project root:  sh 06_resubmission_20261008/code/supplement/run_all.sh
set -e
PY=_workspace/envs/celltypist/bin/python
D=06_resubmission_20261008/code/supplement
$PY -I $D/build_S1_S13_tables.py          # S1-S13 Table parts (10-23 layout)
$PY -I $D/update_S4_lung_family.py        # 10-28: lung disease-scope BH family without CHP
$PY -I $D/update_S5_mKH.py                # 10-28: modified Hartung-Knapp row
$PY -I $D/update_parts_20261027_28.py     # 10-27/28: GPL6480, lung QC sheets, MAGMA windows, blood table, spelling
$PY -I $D/append_S6g_alveolar_macrophage.py  # 10-28: alveolar macrophage tests
$PY -I $D/update_S6_region_label.py         # 10-29: SEA-AD region shown as DLPFC
$PY -I $D/update_lung_donor_merge_20261029.py  # 10-29: GSE302339 cood8 = copd8, lung rows of S6, S12, S13 for 55 donors
$PY -I $D/update_parts_20261030.py        # 10-30: GSE48350 sigQC without GSM1176257 (S3), S4 (c) and S12 (a) legends
$PY -I $D/update_parts_20261105.py        # 11-05: Table D (e), secondary comparisons within GSE47460
$PY -I $D/update_figure_numbers_20261105.py  # 11-05: main-figure citations after renumbering to Fig 1-8
$PY -I $D/trim_parts_20261112.py           # 11-12: keep only tables and rows that reported results depend on
LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript $D/build_S1_S2_Fig.R
$PY -I $D/build_S1_File_methods.py        # 10-29: Supplementary Methods pages of S1 File
$PY -I $D/build_S1_S2_File.py             # merge into S1 File (PDF) and S2 File (xlsx)
$PY -I $D/build_S1_Data.py                # S1 Data
$PY -I $D/set_supplement_fonts.py           # 11-09: Times New Roman in S2 File and S1 Data, as in the main text

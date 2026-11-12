# COPD lung and Alzheimer's disease brain co-expression analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![R 4.5.3](https://img.shields.io/badge/R-4.5.3-276DC3)
![Python 3.9](https://img.shields.io/badge/Python-3.9-3776AB)

Code for **Co-expression network analysis of chronic obstructive pulmonary disease lung and Alzheimer's disease brain reveals a shared myeloid gene program enriched for Alzheimer's disease genetic risk** (Feijie Yao, PLOS ONE submission). All inputs are public; the Axis A and Axis B gene lists are in [`resources/`](resources/).

## Find the code

| By analysis | By figure |
| --- | --- |
| [`scripts/`](scripts/): core pipeline<br>[`analyses/`](analyses/README.md): further analyses | [Fig 1](by_figure/Fig1.md) · [Fig 2](by_figure/Fig2.md) · [Fig 3](by_figure/Fig3.md) · [Fig 4](by_figure/Fig4.md) · [Fig 5](by_figure/Fig5.md) · [Fig 6](by_figure/Fig6.md) · [Fig 7](by_figure/Fig7.md) · [Fig 8](by_figure/Fig8.md) · [Supporting information](by_figure/Supporting_information.md) |

## Repository layout

```text
scripts/        core pipeline by analysis step (run with run_all.sh)
R/  config/     shared functions, paths and stages of the core pipeline
validation/     reference outputs and comparison script
analyses/       further analyses, one folder each
manuscript/     figure code (main_figures/, figure_sources/) and supporting information (supplement/)
by_figure/      figure panels linked to code
resources/      Axis A/B gene lists, NCBI Gene snapshot
env/  data/     environments; input manifest
```

### Core pipeline

| Folder | Analysis | Used in |
| --- | --- | --- |
| [`00_download/`](scripts/00_download/) | Inputs, annotations, Hallmark sets, GWAS preparation | – |
| [`01_wgcna/`](scripts/01_wgcna/) | WGCNA networks, module–trait correlations, module preservation | Fig 2, Fig 3a, Fig A |
| [`02_module_overlap/`](scripts/02_module_overlap/) | Module overlap, Axis A and Axis B | Fig 2, Table A–B |
| [`03_enrichment/`](scripts/03_enrichment/) | GO, KEGG and Hallmark enrichment | Fig 6a, Table H, Fig B |
| [`04_genetics_blood/`](scripts/04_genetics_blood/) | Peripheral blood, MAGMA, disease scope | Fig 6c, Fig 7, Table D, J, K |
| [`05_single_cell/`](scripts/05_single_cell/) | Single-cell QC, annotation, UCell localization | Fig 4, Fig 5b, Table F |
| [`06_virtual_knockout/`](scripts/06_virtual_knockout/) | CellOracle virtual knockout | Fig 8, Table L |

Lung single-cell donor-level results in the article use 55 donors ([`analyses/single_cell_lung_55_donors`](analyses/single_cell_lung_55_donors/)); the core pipeline uses the 56 subject IDs of GEO.

## Running

```sh
python3 scripts/00_download/01_fetch.py   # inputs listed in data/manifest.tsv
sh run_all.sh wgcna                         # one stage at a time, see config/stages.json
sh run_all.sh validate                      # compare with validation/reference/
```

Scripts in `analyses/` and `manuscript/` are run from the project root ([paths](analyses/README.md#paths-in-the-scripts)).

## Environment

R 4.5.3 ([`env/renv.lock`](env/renv.lock)) and Python 3.9.6 ([`env/`](env/)). Key versions: WGCNA 1.74, limma 3.66.0, edgeR 4.8.2, singscore 1.30.0, metafor 5.0.1, UCell 2.14.0, sigQC 0.1.24, decoupleR 2.16.0, CellOracle 0.20.0, MAGMA v1.10 (<https://cncr.nl/research/magma/>). CIBERSORTx was run through its web tool (<https://cibersortx.stanford.edu/>).

## Seeds

WGCNA 20260624; preservation 20260819 (lung), 20260712 (brain); single-cell random gene sets 20260901 (lung), 20260831 (brain); sigQC 20261014; matched random gene sets 20261019; CellOracle 123; UMAP 0.

## Citation and license

Please cite the article ([`CITATION.cff`](CITATION.cff)). Code: MIT ([`LICENSE`](LICENSE)); input data remain under the terms of their original repositories.

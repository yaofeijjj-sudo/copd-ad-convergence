# R7 conditional MAGMA: Axis A AD enrichment conditional on a general microglia gene set (design, frozen 2026-10-16 before any run)

Authorised by the author 2026-10-16 ("这两点现在推进"). Written before any conditional result was seen.

## Question

Is the AD MAGMA enrichment of Axis A (137 genes) explained by Axis A being a microglia/myeloid gene set? Alternative R7 (`paper/研究逻辑_20261009.md` §8): "myeloid genes are AD risk genes anyway".

## Estimand, unit, test

- Estimand: Axis A MAGMA competitive gene-set coefficient (one-sided, beta > 0) in the AD GWAS, with membership in an external general microglia gene set included as a covariate.
- Unit: gene (18,224 genes in the AD gene analysis).
- Software: MAGMA v1.10, `--model condition=MicrogliaRef` (manual v1.10 "Conditional and joint analysis"; de Leeuw et al. 2018 Nat Commun; same mechanism as Jansen et al. 2019 Nat Genet conditioning on an APOE gene set). Default internal corrections unchanged.
- Gene-level input: frozen SYN01 `GWAS和外周血/SYN01_Entrez_downstream_20261001/results/magma/AD_total_N.genes.raw` (GCST90704648, no biobank, total N). No new gene analysis.
- Gene sets: SYN01 `locked_axes.sets` (Axis A 137, Axis B 434) plus the reference set below, in one set-annot file.

## Reference set (fixed before running)

- Planned: Patir et al. 2019 Glia 249-gene human microglia core signature. The supplementary file could not be retrieved (Wiley Cloudflare check; not bypassed). Author-approved fallback used instead.
- Used: BRETIGEA 1.0.3 `markers_df_human_brain`, cell = "mic" (McKenzie et al. 2018 Sci Rep; human consensus microglia markers ranked by specificity; external to this project).
  - Primary: top 250 markers (size chosen to match the planned 249-gene reference, fixed now).
  - Sensitivity: all 1,000 microglia markers.
- Symbol -> Entrez: org.Hs.eg.db 3.22.0, official symbols only (no aliases), one-to-many dropped; same rule as the deconvolution branch.
- Axis A genes are NOT removed from the reference set; overlap is reported. MAGMA conditional regression distributes shared association between the two indicators (manual); this is the standard use.

## Models

1. Unconditional reference: MicrogliaRef alone (descriptive: does a general microglia set enrich AD?).
2. Primary: AxisA and AxisB each conditional on MicrogliaRef (top 250).
3. Sensitivity: same with the 1,000-marker reference.
4. Descriptive: Axis A members outside the reference set ("AxisA_nonRef"), unconditional.

Multiplicity: Holm across the two axes within each model family (as in the APOE sensitivity). No selection between reference sizes; both reported.

## Interpretation rules (written before results)

- Axis A conditional Holm P < 0.05: Axis A's AD enrichment is not fully explained by general microglia membership.
- Otherwise: Axis A's AD enrichment cannot be distinguished from a general microglia background; the manuscript states only that Axis A is a microglia/myeloid response associated with AD genetic risk (the current wording already does not claim specificity).
- If primary and sensitivity disagree, report both; primary governs the manuscript sentence.

## Risks

- BRETIGEA markers were derived from expression specificity, not genetics; top-N cut is arbitrary but fixed here.
- Conditioning on a set that shares many Axis A genes reduces power; a null result may reflect overlap rather than absence of extra signal (report the overlap).
- Association-level only; no causal inference.

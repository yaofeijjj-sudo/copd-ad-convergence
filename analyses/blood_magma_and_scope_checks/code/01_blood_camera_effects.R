# Advice 3.1 (DESIGN.md): CAMERA on the fixed Axis A set, subject-level Axis A score differences, AddNeuroMed duplicate check.
# Data loading and models reuse the SYN01 / 10-16 blood code; baseline t must equal the frozen SYN01 statistics.
# Run from project root:
#   LC_ALL=en_US.UTF-8 Rscript 03_candidate_analysis/AxisA_auxiliary_checks_20261019/code/01_blood_camera_effects.R
.libPaths(c("_workspace/r_libs_blood_deconv_20261016", "_workspace/r_libs", .libPaths()))
suppressPackageStartupMessages({ library(AnnotationDbi); library(Biobase); library(data.table); library(GEOquery); library(limma); library(org.Hs.eg.db) })
options(stringsAsFactors = FALSE)
ROOT <- normalizePath(getwd())
O <- file.path(ROOT, "03_candidate_analysis/AxisA_auxiliary_checks_20261019")
SYN <- file.path(ROOT, "GWAS和外周血/SYN01_Entrez_downstream_20261001")
OLD <- file.path(ROOT, "GWAS和外周血/current_figure_package/six_blood_models_candidate_20260927")
dir.create(file.path(O, "sandbox/GWAS和外周血/SYN01_Entrez_downstream_20261001/results"), recursive = TRUE, showWarnings = FALSE)
src <- parse(file.path(SYN, "code/refit_gse248493.R"))
for (e in src) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) eval(e, globalenv())
raw_dir <- file.path(ROOT, "02_data/raw/bulk_blood_external")
gse248493_counts_file <- file.path(OLD, "inputs/GSE248493_merged_counts_20260707.csv.gz")
project_root <- file.path(O, "sandbox")
axis <- as.character(fread(file.path(ROOT, "03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisA_confirmed_genes.csv"))$entrez_id)
stopifnot(length(axis) == 137)

# each entry: y (matrix or voom EList), design, contrast vector, t (for the gate)
M <- list(); prep_keep <- list()
from_syn <- function(prepared, sample_dt) {
  r <- run_limma_contrast(prepared, sample_dt, "age_gender_adjusted", c("age", "gender"))
  y <- if (identical(prepared$data_type, "rnaseq_counts")) readRDS(file.path(project_root, "GWAS和外周血/SYN01_Entrez_downstream_20261001/results/GSE248493_voom.rds")) else r$expr_sub
  list(y = y, design = r$design, contrast = r$contrast, t = r$gene_stats[, .(entrez_id, t)])
}
p <- prepare_microarray("GSE54837"); M$GSE54837 <- from_syn(p, sample_table_gse54837(p))
p <- prepare_gse248493_counts(gse248493_counts_file); M$GSE248493 <- from_syn(p, sample_table_gse248493(p))
p <- prepare_microarray("GSE63060"); M$GSE63060 <- from_syn(p, sample_table_addneuromed(p, "GSE63060_AD_blood_vs_CTL_casecontrol")); prep_keep$GSE63060 <- p
p <- prepare_microarray("GSE63061"); M$GSE63061 <- from_syn(p, sample_table_addneuromed(p, "GSE63061_AD_blood_vs_CTL_casecontrol")); prep_keep$GSE63061 <- p
rm(p); gc()
local({  # E-MTAB-5278, as refit_emtab5278.R / 10-16 code
  raw <- file.path(raw_dir, "E-MTAB-5278_20260929")
  s <- fread(file.path(ROOT, "GWAS和外周血/current_figure_package/replacement_age_sex_EMTAB5278_20260929/sample_manifest_emtab5278.csv"))[include_primary == TRUE]; stopifnot(nrow(s) == 109L)
  s[, group := factor(group, levels = c("Control", "Case"))]; s[, sex := factor(sex, levels = c("F", "M"))]
  mf <- file.path(raw, "H1_QASMC_Blood_Hs_MRNA_normMatrix.txt"); header <- strsplit(readLines(mf, n = 1L), "\t", fixed = TRUE)[[1]]
  m <- fread(mf, skip = 2L, header = FALSE); setnames(m, header); setnames(m, "Hybridization Name", "gene_symbol")
  symbol <- toupper(trimws(as.character(m$gene_symbol))); symbol[!nzchar(symbol) | symbol %in% c("NA", "N/A", "---")] <- NA_character_
  se <- as.matrix(m[, s$array_file, with = FALSE]); storage.mode(se) <- "double"
  idl <- fread(file.path(SYN, "results/EMTAB5278_identity_ledger.csv"), colClasses = "character")
  entrez <- idl$entrez_id[match(symbol, idl$input_symbol)]; valid <- !is.na(entrez) & nzchar(entrez)
  expr <- limma::avereps(se[valid, , drop = FALSE], ID = entrez[valid])
  d <- model.matrix(~ 0 + group + age + sex, data = s); rownames(d) <- s$array_file
  ct <- limma::makeContrasts(groupCase - groupControl, levels = d)[, 1]
  f <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr, d), ct))
  M[["E-MTAB-5278"]] <<- list(y = expr, design = d, contrast = ct, t = as.data.table(limma::topTable(f, 1, Inf, sort.by = "none"), keep.rownames = "entrez_id")[, .(entrez_id, t)])
})
local({  # GSE248417, as refit_gse248417.R / 10-16 code
  raw <- file.path(raw_dir, "GSE248417_20260927")
  manifest <- fread(file.path(OLD, "sample_manifest_gse248417.csv")); manifest[, include_primary := tolower(as.character(include_primary)) == "true"]
  s <- manifest[include_primary == TRUE]; stopifnot(nrow(s) == 97L)
  s[, group := factor(group, levels = c("Control", "AD"))]; s[, sex := factor(sex, levels = c("F", "M"))]
  src <- fread(file.path(raw, "GSE248417_ADvsControl.mRNA.gene_level.deg.txt.gz")); src[, ensembl_id := sub("\\..*$", "", gene_id)]
  mp <- fread(file.path(SYN, "results/GSE248417_Ensembl_Entrez_ledger.csv"), colClasses = "character")[n_ids == "1"]
  src[, entrez_id := mp$ENTREZID[match(ensembl_id, mp$ENSEMBL)]]; cols <- s$sample_id
  cd <- src[!is.na(entrez_id), lapply(.SD, sum), by = entrez_id, .SDcols = cols]
  counts <- as.matrix(cd[, ..cols]); storage.mode(counts) <- "double"; rownames(counts) <- cd$entrez_id
  d <- model.matrix(~ 0 + group + age + sex, data = s); rownames(d) <- s$sample_id
  ct <- limma::makeContrasts(groupAD - groupControl, levels = d)[, 1]
  cpm <- t(t(counts) / colSums(counts) * 1e6); keep <- rowSums(cpm >= 1) >= max(3L, floor(min(table(s$group)) / 2))
  v <- limma::voom(edgeR::calcNormFactors(edgeR::DGEList(counts[keep, ]), method = "TMM"), d, plot = FALSE)
  f <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v, d), ct))
  M$GSE248417 <<- list(y = v, design = d, contrast = ct, t = as.data.table(limma::topTable(f, 1, Inf, sort.by = "none"), keep.rownames = "entrez_id")[, .(entrez_id, t)])
})
ids <- c("GSE54837", "GSE248493", "E-MTAB-5278", "GSE63060", "GSE63061", "GSE248417"); M <- M[ids]

# ---- gate ----
gate <- rbindlist(lapply(ids, function(id) {
  fz <- fread(file.path(SYN, "results", paste0("ranked_gene_stats_", id, ".csv.gz")), colClasses = c(entrez_id = "character"))
  b <- merge(M[[id]]$t, fz[, .(entrez_id, t_frozen = t)], by = "entrez_id")
  data.table(dataset_id = id, n = nrow(M[[id]]$t), n_frozen = nrow(fz), n_matched = nrow(b), max_abs_t_diff = max(abs(b$t - b$t_frozen))) }))
print(gate); stopifnot(all(gate$n == gate$n_frozen), all(gate$n_matched == gate$n_frozen), all(gate$max_abs_t_diff < 1e-8))
fwrite(gate, file.path(O, "qa/blood_baseline_gate.tsv"), sep = "\t")

# ---- CAMERA + subject-level score ----
meta <- data.table(dataset_id = ids, disease = c("COPD", "COPD", "COPD", "AD", "AD", "AD"), sample_type = c("whole blood", "PBMC", "whole blood", "whole blood", "whole blood", "whole blood"))
fgs <- fread(file.path(SYN, "results/fgsea_six_results.csv"))
res <- rbindlist(lapply(ids, function(id) {
  m <- M[[id]]; genes <- rownames(if (is(m$y, "EList")) m$y$E else m$y); idx <- which(genes %in% axis)
  cm <- limma::camera(m$y, index = idx, design = m$design, contrast = m$contrast, inter.gene.cor = NA)
  E <- if (is(m$y, "EList")) m$y$E else m$y
  sc <- singscore::simpleScore(singscore::rankGenes(E), upSet = genes[idx])$TotalScore; z <- as.numeric(scale(sc))
  f <- lm(z ~ 0 + m$design); b <- coef(f); V <- vcov(f); k <- m$contrast
  est <- sum(k * b); se <- sqrt(drop(t(k) %*% V %*% k)); df <- f$df.residual
  data.table(dataset_id = id, n = nrow(m$design), axis_present = length(idx), camera_direction = cm$Direction, camera_P = cm$PValue, camera_corr = cm$Correlation,
             score_diff_SD = est, ci_low = est - qt(.975, df) * se, ci_high = est + qt(.975, df) * se, score_P = 2 * pt(-abs(est / se), df)) }))
res[, camera_BH := p.adjust(camera_P, "BH")]
res <- merge(meta, res, by = "dataset_id"); res <- merge(res, fgs[, .(dataset_id, fgsea_NES, fgsea_p)], by = "dataset_id", all.x = TRUE)
res <- res[match(ids, dataset_id)]
fwrite(res, file.path(O, "results/blood_camera_and_effects.tsv"), sep = "\t"); print(res)

# ---- AddNeuroMed cross-dataset duplicate check ----
a <- prep_keep$GSE63060$expr_gene; b <- prep_keep$GSE63061$expr_gene; g <- intersect(rownames(a), rownames(b))
ac <- a[g, ] - rowMeans(a[g, ]); bc <- b[g, ] - rowMeans(b[g, ]); R <- cor(ac, bc)
best_b <- apply(R, 1, function(x) { o <- order(x, decreasing = TRUE); c(colnames(R)[o[1]], x[o[1]], x[o[2]]) })
du <- data.table(GSE63060 = rownames(R), GSE63061 = best_b[1, ], r_best = as.numeric(best_b[2, ]), r_second = as.numeric(best_b[3, ]))
mutual <- apply(R, 2, which.max); du[, mutual := mutual[GSE63061] == match(GSE63060, rownames(R))]
du[, flagged := mutual & (r_best - r_second) >= 0.1]
cat("AddNeuroMed: genes", length(g), "| samples", nrow(R), "x", ncol(R), "| median best r", round(median(du$r_best), 3), "| flagged", sum(du$flagged), "\n")
fwrite(du[order(-r_best)], file.path(O, "results/addneuromed_cross_dataset_matches.tsv"), sep = "\t")
print(head(du[order(-(r_best - r_second))], 5))
writeLines(capture.output(sessionInfo()), file.path(O, "qa/sessionInfo_01.txt"))

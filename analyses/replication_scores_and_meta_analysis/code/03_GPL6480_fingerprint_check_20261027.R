# GPL6480 donors vs lung discovery cohort GSE76925: metadata fingerprints and expression similarity (2026-10-27). Run from project root.
suppressPackageStartupMessages({ library(data.table) })
ROOT <- normalizePath("."); P <- function(x) file.path(ROOT, x)
# helpers copied from AxisA_singscore_official_rerun_20261004/code/02_run_specificity.R
hist <- fread(cmd = "gzip -cd 03_candidate_analysis/WGCNA_Entrez_training_only_20260930/references/gene_history.tsv.gz", colClasses = "character")
setnames(hist, 1, "tax_id"); hist <- hist[tax_id == "9606" & GeneID != "-"]; repl <- setNames(hist$GeneID, hist$Discontinued_GeneID)
current_id <- function(x) { for (i in 1:5) { hit <- !is.na(x) & x %in% names(repl); if (!any(hit)) break; x[hit] <- repl[x[hit]] }; x }
read_matrix <- function(file) {
  tab <- fread(cmd = paste("gzip -cd", shQuote(P(file)), "| sed -n '/^!series_matrix_table_begin/,$p' | sed '1d;/^!series_matrix_table_end/,$d'"), check.names = FALSE)
  m <- as.matrix(tab[, -1L, with = FALSE]); storage.mode(m) <- "numeric"; rownames(m) <- as.character(tab[[1L]]); m
}
keep_max_mad <- function(m, gene) {
  ok <- !is.na(gene); m <- m[ok, , drop = FALSE]; gene <- gene[ok]
  v <- apply(m, 1, mad, na.rm = TRUE); o <- order(gene, -v); o <- o[is.finite(v[o])]
  o <- o[!duplicated(gene[o])]; out <- m[o, , drop = FALSE]; rownames(out) <- gene[o]
  out[rowSums(is.finite(out)) == ncol(out), , drop = FALSE]
}
a <- fread("02_data/raw/platform_annotations/GPL6480_annot_table.csv", colClasses = "character", select = c("ID", "Gene ID")); setnames(a, c("probe", "gene"))
m <- read_matrix("02_data/raw/bulk_lung_copd/GSE47460_GPL6480/GSE47460-GPL6480_series_matrix.txt.gz")
a <- a[probe %in% rownames(m) & !is.na(gene) & nzchar(gene) & !grepl("/", gene)]
m <- keep_max_mad(m[a$probe, , drop = FALSE], current_id(a$gene))
ph <- fread("GWAS_and_blood/analysis/05_disease_specificity_gwas/results/axisA_disease_context_selectivity_20260711/metadata/GSE47460_GPL6480_pheno_raw.csv")
ph <- ph[match(colnames(m), geo_accession)]; ph[, donor := sub("^(LT[0-9]{6}).*$", "\\1", title)]
ph <- ph[`disease state:ch1` %in% c("Control", "Chronic Obstructive Lung Disease")]
ids <- unique(ph$donor); b <- sapply(ids, function(id) rowMeans(m[, ph[donor == id]$geo_accession, drop = FALSE])); rownames(b) <- rownames(m)
coh <- readRDS("03_candidate_analysis/AxisA_member_direction_robustness_20261018/source_data/cohorts.rds")$cohorts
e <- coh$GSE76925$expr[, coh$GSE76925$dat$sample]
g <- intersect(rownames(e), rownames(b)); za <- t(scale(t(e[g, ]))); zb <- t(scale(t(b[g, ]))); ok <- rowSums(!is.finite(za)) == 0 & rowSums(!is.finite(zb)) == 0
cc <- cor(za[ok, ], zb[ok, ]); cat("genes", sum(ok), "| GSE76925", nrow(cc), "| GPL6480 donors", ncol(cc), "\n")
# fingerprint pairs (age, sex, FEV1 % predicted pre- or post-bronchodilator)
e1 <- new.env(); load("WGCNA_figures/reproduction_materials/analysis/01_lung_wgcna/results/COPD_WGCNA/copd_wgcna_final_env.RData", envir = e1)
ct <- as.data.table(e1$clinical_train); stopifnot(nrow(ct) == nrow(cc))
d76 <- data.table(s76 = rownames(cc), age = round(as.numeric(ct$age)), sex = fifelse(as.character(ct$sex) %in% c("1", "M", "Male"), "M", "F"), fev1 = round(as.numeric(ct$FEV1_pp)))
stopifnot(identical(d76$s76, as.character(coh$GSE76925$dat$sample)))
d64 <- ph[!duplicated(donor), .(donor, age = round(as.numeric(`age:ch1`)), sex = fifelse(grepl("Male", `Sex:ch1`) & !grepl("Female", `Sex:ch1`), "M", "F"),
          pre = round(as.numeric(`%predicted fev1 (pre-bd):ch1`)), post = round(as.numeric(`%predicted fev1 (post-bd):ch1`)))]
fp <- unique(rbind(merge(d64, d76, by.x = c("age", "sex", "pre"), by.y = c("age", "sex", "fev1"))[, .(s76, donor, age, sex, fev1 = pre)],
                   merge(d64, d76, by.x = c("age", "sex", "post"), by.y = c("age", "sex", "fev1"))[, .(s76, donor, age, sex, fev1 = post)]))
fp[, `:=`(cor = cc[cbind(s76, donor)], rank = mapply(function(s, d) rank(-cc[s, ])[d], s76, donor), n = ncol(cc),
          best_other = mapply(function(s, d) max(cc[s, colnames(cc) != d]), s76, donor))]
print(fp); cat("median best r of GSE76925 samples:", round(median(apply(cc, 1, max)), 3), "\n")
fwrite(fp, "03_candidate_analysis/AxisA_singscore_official_rerun_20261004/qa/GPL6480_vs_GSE76925_fingerprint_pairs_20261027.tsv", sep = "\t")

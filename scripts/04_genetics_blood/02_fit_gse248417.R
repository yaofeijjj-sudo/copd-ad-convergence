source("R/runtime.R")


args <- commandArgs(trailingOnly = TRUE)
root_arg <- grep("^--root=", args, value = TRUE)
project_root <- if (length(root_arg)) sub("^--root=", "", root_arg[[1]]) else getwd()
project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
.libPaths(c(file.path(project_root, "env/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(fgsea)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
set.seed(params$blood_rnaseq_preprocessing_seed)

old <- file.path(project_root, "results/blood_ranks")
out <- file.path(project_root, "results/genetics_blood/results")
raw <- file.path(project_root, "data/raw/bulk_blood_external/GSE248417")
matrix_file <- file.path(raw, "GSE248417_ADvsControl.mRNA.gene_level.deg.txt.gz")

axis_entrez <- as.character(fread(file.path(project_root, "results/module_overlap/results/cross/AxisA_confirmed_genes.csv"))$entrez_id)
axis_symbols <- axis_entrez
manifest <- fread(file.path(project_root, "results/blood_metadata/sample_manifest_gse248417.csv"))
stopifnot(nrow(manifest) == 98L, uniqueN(manifest$sample_id) == 98L)
manifest[, include_primary := tolower(as.character(include_primary)) == "true"]
manifest[, group := factor(group, levels = c("Control", "AD"))]
manifest[, sex := factor(sex, levels = c("F", "M"))]
stopifnot(sum(manifest$group == "AD") == 49L, sum(manifest$group == "Control") == 49L)
stopifnot(sum(manifest$include_primary) == 97L)

source <- fread(matrix_file)
stopifnot(nrow(source) == 18891L, uniqueN(source$gene_id) == nrow(source))
sample_cols <- manifest$sample_id
stopifnot(all(sample_cols %in% names(source)))
source[, symbol := gene_name]
source[, ensembl_id := sub("\\..*$", "", gene_id)]
keys <- intersect(unique(source$ensembl_id), AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENSEMBL"))
mp <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = keys, keytype = "ENSEMBL", columns = "ENTREZID"))
mp <- unique(mp[!is.na(ENTREZID)])
mp[, n_ids := uniqueN(ENTREZID), by = ENSEMBL]
fwrite(mp, file.path(out, "GSE248417_Ensembl_Entrez_ledger.csv"))
source[, entrez_id := as.character(mp[n_ids == 1]$ENTREZID[match(ensembl_id, mp[n_ids == 1]$ENSEMBL)])]
mapping_qc <- data.table(
  input_gene_rows = nrow(source),
  unique_symbols = uniqueN(source$symbol[!is.na(source$symbol)]),
  mapped_rows = sum(!is.na(source$entrez_id)),
  mapped_entrez = uniqueN(source$entrez_id[!is.na(source$entrez_id)]),
  axis_input_symbols = length(axis_symbols),
  axis_mapped_entrez = length(axis_entrez)
)

usable <- source[!is.na(entrez_id)]
counts_dt <- usable[, lapply(.SD, sum), by = entrez_id, .SDcols = sample_cols]
counts <- as.matrix(counts_dt[, ..sample_cols])
storage.mode(counts) <- "double"
rownames(counts) <- counts_dt$entrez_id
stopifnot(!anyNA(counts), all(is.finite(counts)), all(counts >= 0), all(colSums(counts) > 0))
fwrite(mapping_qc, file.path(out, "gene_mapping_qc.csv"))

run_model <- function(name, selected, covariates) {
  samples <- manifest[selected]
  stopifnot(!anyNA(samples$group))
  if (length(covariates)) {
    stopifnot(all(vapply(samples[, ..covariates], function(x) !anyNA(x), logical(1))))
  }
  counts_sub <- counts[, samples$sample_id, drop = FALSE]
  samples[, group := factor(group, levels = c("Control", "AD"))]
  formula <- as.formula(paste("~ 0 +", paste(c("group", covariates), collapse = " + ")))
  design <- model.matrix(formula, data = samples)
  rownames(design) <- samples$sample_id
  stopifnot(qr(design)$rank == ncol(design))
  contrast <- limma::makeContrasts(groupAD - groupControl, levels = design)

  cpm <- t(t(counts_sub) / colSums(counts_sub) * 1e6)
  threshold_n <- max(3L, floor(min(table(samples$group)) / 2))
  keep <- rowSums(cpm >= 1) >= threshold_n
  stopifnot(sum(keep) > 1000L)
  filtered <- counts_sub[keep, , drop = FALSE]
  dge <- edgeR::calcNormFactors(edgeR::DGEList(filtered), method = "TMM")
  fwrite(as.data.table(dge$samples, keep.rownames = "sample_id"), file.path(out, "GSE248417_TMM_factors.csv"))
  v <- limma::voom(dge, design, plot = FALSE, save.plot = TRUE)
  saveRDS(v, file.path(out, "GSE248417_voom.rds"))
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v, design), contrast))
  stats <- as.data.table(limma::topTable(fit, coef = 1, number = Inf, sort.by = "none"), keep.rownames = "entrez_id")
  stopifnot(!anyNA(stats$t), uniqueN(stats$entrez_id) == nrow(stats))
  setorder(stats, -t)
  fwrite(stats, file.path(out, paste0("gene_stats_", name, ".csv.gz")))

  data.table(n_AD = sum(samples$group == "AD"), n_Control = sum(samples$group == "Control"), gene_rows_after_filter = nrow(filtered))
}

primary <- run_model("primary_age_sex_97", manifest$include_primary, c("age", "sex"))
fwrite(primary, file.path(out, "GSE248417_model_qc.csv"))
capture.output(sessionInfo(), file = file.path(out, "session_info_GSE248417_refit.txt"))
print(primary)

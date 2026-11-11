source("R/runtime.R")
args <- commandArgs(trailingOnly = TRUE)
root_arg <- grep("^--root=", args, value = TRUE)
root <- if (length(root_arg)) sub("^--root=", "", root_arg[[1]]) else getwd()
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
.libPaths(c(file.path(root, "env/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(limma)
  library(fgsea)
})
stopifnot(
  as.character(packageVersion("limma")) == "3.66.0",
  as.character(packageVersion("fgsea")) == "1.36.2"
)

out <- file.path(root, "results/genetics_blood/results")
old <- file.path(root, "results/blood_ranks")
raw <- file.path(root, "data/raw/bulk_blood_external/E-MTAB-5278")
matrix_file <- file.path(raw, "H1_QASMC_Blood_Hs_MRNA_normMatrix.txt")
manifest_file <- file.path(root, "results/blood_metadata/sample_manifest_emtab5278.csv")

axis_entrez <- as.character(fread(file.path(root, "results/module_overlap/results/cross/AxisA_confirmed_genes.csv"))$entrez_id)
axis_symbols <- axis_entrez
manifest <- fread(manifest_file)
stopifnot(
  nrow(manifest) == 234L,
  sum(manifest$include_primary) == 109L
)
samples <- manifest[include_primary == TRUE]
stopifnot(
  nrow(samples) == 109L,
  uniqueN(samples$donor_id) == 109L,
  uniqueN(samples$array_file) == 109L,
  sum(samples$group == "Case") == 51L,
  sum(samples$group == "Control") == 58L,
  all(samples$smoking_status == "smoker"),
  !anyNA(samples[, .(group, age, sex)])
)
samples[, group := factor(group, levels = c("Control", "Case"))]
samples[, sex := factor(sex, levels = c("F", "M"))]
design <- model.matrix(~ 0 + group + age + sex, data = samples)
rownames(design) <- samples$array_file
stopifnot(qr(design)$rank == ncol(design), ncol(design) == 4L)
contrast <- limma::makeContrasts(groupCase - groupControl, levels = design)

header <- strsplit(readLines(matrix_file, n = 1L), "\t", fixed = TRUE)[[1]]
stopifnot(
  header[[1]] == "Hybridization Name",
  length(header) == 225L, uniqueN(header) == 225L,
  all(samples$array_file %in% header)
)
matrix_dt <- fread(matrix_file, skip = 2L, header = FALSE)
stopifnot(ncol(matrix_dt) == length(header))
setnames(matrix_dt, header)
setnames(matrix_dt, "Hybridization Name", "gene_symbol")
symbol <- toupper(trimws(as.character(matrix_dt$gene_symbol)))
symbol[!nzchar(symbol) | symbol %in% c("NA", "N/A", "---")] <- NA_character_
sample_cols <- samples$array_file
source_expr <- as.matrix(matrix_dt[, ..sample_cols])
storage.mode(source_expr) <- "double"
stopifnot(
  nrow(source_expr) > 10000L,
  !anyNA(source_expr), all(is.finite(source_expr)),
  max(source_expr) < 100, min(source_expr) > -100
)
colnames(source_expr) <- samples$array_file
stopifnot(identical(colnames(source_expr), rownames(design)))

keys <- unique(na.omit(symbol))
mp <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = intersect(keys, AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")), keytype = "SYMBOL", columns = "ENTREZID"))
mp <- unique(mp[!is.na(ENTREZID)])
mp[, n_ids := uniqueN(ENTREZID), by = SYMBOL]
canon <- mp[n_ids == 1, .(input_symbol = SYMBOL, entrez_id = as.character(ENTREZID), mapping_rule = "unique_SYMBOL")]
missing <- setdiff(keys, canon$input_symbol)
ap <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = intersect(missing, AnnotationDbi::keys(org.Hs.eg.db, keytype = "ALIAS")), keytype = "ALIAS", columns = "ENTREZID"))
ap <- unique(ap[!is.na(ENTREZID)])
ap[, n_ids := uniqueN(ENTREZID), by = ALIAS]
identities <- rbind(canon, ap[n_ids == 1, .(input_symbol = ALIAS, entrez_id = as.character(ENTREZID), mapping_rule = "unique_official_historical_symbol")])
fwrite(identities, file.path(out, "EMTAB5278_identity_ledger.csv"))
entrez <- identities$entrez_id[match(symbol, identities$input_symbol)]
valid <- !is.na(entrez) & nzchar(entrez)
expr <- limma::avereps(source_expr[valid, , drop = FALSE], ID = entrez[valid])
storage.mode(expr) <- "double"
stopifnot(!anyNA(expr), all(is.finite(expr)), nrow(expr) > 10000L)
covered <- intersect(axis_entrez, rownames(expr))
mapping_qc <- data.table(
  matrix_feature_rows = nrow(matrix_dt),
  matrix_sample_columns = ncol(matrix_dt) - 1L,
  selected_cases = sum(samples$group == "Case"),
  selected_controls = sum(samples$group == "Control"),
  unique_selected_donors = uniqueN(samples$donor_id),
  mapped_feature_rows = sum(valid),
  mapped_entrez_genes = nrow(expr),
  axis_input_symbols = length(axis_symbols),
  axis_mapped_entrez = length(axis_entrez),
  axis_covered = length(covered),
  design_rank = qr(design)$rank,
  design_columns = ncol(design),
  min_expression = min(source_expr),
  max_expression = max(source_expr)
)
fwrite(mapping_qc, file.path(out, "gene_mapping_qc_emtab5278.csv"))
stopifnot(length(covered) >= 119L)

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr, design), contrast))
stats <- as.data.table(limma::topTable(fit,
  coef = 1L,
  number = Inf, sort.by = "none"
), keep.rownames = "entrez_id")
stopifnot(
  uniqueN(stats$entrez_id) == nrow(stats),
  !anyNA(stats$t), all(is.finite(stats$t))
)
stats[, `:=`(
  dataset_id = "E-MTAB-5278", model_id = "age_sex_adjusted",
  rank_metric = "limma moderated t; positive = COPD higher than current-smoker control"
)]
setorder(stats, -t)
fwrite(stats, file.path(out, "gene_stats_emtab5278.csv.gz"))

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_EMTAB5278.txt"))

source("R/runtime.R")
args <- commandArgs(trailingOnly = TRUE)
root_arg <- grep("^--root=", args, value = TRUE)
project_root <- if (length(root_arg)) sub("^--root=", "", root_arg[[1]]) else getwd()
project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
local_r_lib <- file.path(project_root, "env/library")
if (dir.exists(local_r_lib)) .libPaths(c(local_r_lib, .libPaths()))
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(Biobase)
  library(data.table)
  library(fgsea)
  library(GEOquery)
  library(limma)
  library(org.Hs.eg.db)
})
options(stringsAsFactors = FALSE)
set.seed(params$blood_preprocessing_seed)
raw_dir <- file.path(project_root, "data/raw/bulk_blood_external")
out_dir <- file.path(project_root, "results/blood_ranks")
gse248493_counts_file <- file.path(out_dir, "inputs/GSE248493_merged_counts.csv.gz")
message_ts <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

normalize_na <- function(x) {
  y <- trimws(as.character(x))
  y[y %in% c("", "NA", "na", "N/A", "n/a", "NaN", "NULL", "null", "Unknown", "unknown")] <- NA_character_
  y
}

first_token <- function(x) {
  x <- normalize_na(x)
  vapply(strsplit(x, "\\s*(///|//|;|,|\\|)\\s*"), function(y) {
    y <- trimws(y)
    y <- y[nzchar(y)]
    if (length(y) == 0) NA_character_ else y[[1]]
  }, character(1))
}

clean_symbol <- function(x) {
  y <- toupper(first_token(x))
  y[y %in% c("NA", "N/A", "NA.", "NAN", "---", "0")] <- NA_character_
  y
}

clean_entrez <- function(x) {
  y <- first_token(x)
  y <- sub("_at$", "", y)
  y <- sub("_s_at$", "", y)
  y <- sub("_x_at$", "", y)
  y[!grepl("^[0-9]+$", y)] <- NA_character_
  y
}

as_numeric_clean <- function(x) {
  suppressWarnings(as.numeric(normalize_na(x)))
}

extract_characteristic <- function(pheno, key) {
  direct_col <- paste0(key, ":ch1")
  if (direct_col %in% colnames(pheno)) {
    return(normalize_na(pheno[[direct_col]]))
  }

  char_cols <- grep("^characteristics", colnames(pheno), value = TRUE)
  out <- rep(NA_character_, nrow(pheno))
  for (col in char_cols) {
    val <- as.character(pheno[[col]])
    hit <- grepl(paste0("^\\s*", key, "\\s*:"), val, ignore.case = TRUE)
    out[hit] <- trimws(sub("^\\s*[^:]+:\\s*", "", val[hit]))
  }
  normalize_na(out)
}

read_gene_set <- function(path, expected_n, label) {
  if (!file.exists(path)) {
    stop("Missing gene set file: ", path)
  }
  dt <- fread(path)
  gene_col <- intersect(c("gene_symbol", "symbol", "Gene", "gene"), names(dt))
  if (!length(gene_col)) {
    stop("No gene column found in ", path)
  }
  genes <- sort(unique(clean_symbol(dt[[gene_col[[1]]]])))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) != expected_n) {
    stop(label, " gene count is ", length(genes), ", expected ", expected_n)
  }
  genes
}

map_symbols_to_entrez <- function(symbols) {
  entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  data.table(
    gene_symbol = names(entrez),
    entrez_id = unname(as.character(entrez)),
    mapped = !is.na(entrez)
  )
}

map_entrez_to_symbols <- function(entrez_ids) {
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = as.character(entrez_ids),
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  out <- unname(as.character(symbols))
  out[is.na(out)] <- as.character(entrez_ids)[is.na(out)]
  out
}

pick_col <- function(cols, candidates) {
  exact <- intersect(candidates, cols)
  if (length(exact)) {
    return(exact[[1]])
  }
  lower_cols <- tolower(cols)
  lower_candidates <- tolower(candidates)
  idx <- match(lower_candidates, lower_cols)
  idx <- idx[!is.na(idx)]
  if (length(idx)) {
    return(cols[[idx[[1]]]])
  }
  NA_character_
}

load_local_eset <- function(gse_id) {
  paths <- c(
    file.path(raw_dir, paste0(gse_id, ""), paste0(gse_id, "_series_matrix.txt.gz")),
    file.path(raw_dir, paste0(gse_id, "_series_matrix.txt.gz"))
  )
  path <- paths[file.exists(paths)][1]
  if (is.na(path)) {
    path <- file.path(raw_dir, paste0(gse_id, "_series_matrix.txt.gz"))
  }
  if (!file.exists(path)) {
    stop("Missing local series matrix for ", gse_id, ": ", path)
  }
  GEOquery::getGEO(filename = path)
}

prepare_microarray <- function(gse_id) {
  message_ts("Loading ", gse_id)
  eset <- load_local_eset(gse_id)
  expr <- Biobase::exprs(eset)
  storage.mode(expr) <- "numeric"
  expression_scale <- if (max(expr, na.rm = TRUE) > 100) {
    expr <- log2(pmax(expr, 1))
    "log2_transformed_by_script"
  } else {
    "already_log2_scale"
  }

  fdata <- Biobase::fData(eset)
  entrez_col <- pick_col(
    colnames(fdata),
    c("ENTREZ_GENE_ID", "Entrez_Gene_ID", "ENTREZID", "Gene ID", "GeneID", "Entrez Gene")
  )
  symbol_col <- pick_col(
    colnames(fdata),
    c("Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol", "SYMBOL", "ILMN_Gene")
  )

  entrez_id <- rep(NA_character_, nrow(fdata))
  if (!is.na(entrez_col)) {
    entrez_id <- clean_entrez(fdata[[entrez_col]])
  }

  symbol <- rep(NA_character_, nrow(fdata))
  if (!is.na(symbol_col)) {
    symbol <- clean_symbol(fdata[[symbol_col]])
  }

  missing_entrez <- is.na(entrez_id) & !is.na(symbol)
  if (any(missing_entrez)) {
    mapped <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = unique(symbol[missing_entrez]),
      keytype = "SYMBOL",
      column = "ENTREZID",
      multiVals = "first"
    )
    entrez_id[missing_entrez] <- unname(as.character(mapped[symbol[missing_entrez]]))
  }

  entrez_id <- entrez_id[match(rownames(expr), rownames(fdata))]
  valid <- !is.na(entrez_id) & nzchar(entrez_id) & rowMeans(is.na(expr)) < 1
  expr_gene <- limma::avereps(expr[valid, , drop = FALSE], ID = entrez_id[valid])
  storage.mode(expr_gene) <- "numeric"
  gene_dt <- data.table(
    entrez_id = rownames(expr_gene),
    gene_symbol = map_entrez_to_symbols(rownames(expr_gene))
  )

  list(
    gse_id = gse_id,
    eset = eset,
    pheno = Biobase::pData(eset),
    expr_gene = expr_gene,
    gene_dt = gene_dt,
    platform_id = Biobase::annotation(eset),
    data_type = "microarray_log_expression",
    expression_scale = expression_scale,
    raw_features = nrow(expr),
    gene_level_universe = nrow(expr_gene),
    entrez_col = entrez_col,
    symbol_col = symbol_col
  )
}

sample_key_from_title <- function(title) {
  key <- sub(".*\\[([^]]+)\\].*", "\\1", as.character(title))
  gsub("\\s+", "_", key)
}

prepare_gse248493_counts <- function(count_file) {
  message_ts("Loading GSE248493 count matrix")
  if (!file.exists(count_file)) {
    stop("Missing merged GSE248493 count matrix: ", count_file)
  }
  eset <- load_local_eset("GSE248493")
  pheno <- Biobase::pData(eset)
  sample_map <- data.table(
    sample_id = rownames(pheno),
    geo_accession = as.character(pheno$geo_accession),
    title = as.character(pheno$title),
    count_sample_id = sample_key_from_title(pheno$title),
    disease = normalize_na(extract_characteristic(pheno, "disease")),
    tissue = normalize_na(extract_characteristic(pheno, "tissue")),
    age = as_numeric_clean(extract_characteristic(pheno, "age")),
    gender = factor(normalize_na(extract_characteristic(pheno, "Sex"))),
    race = normalize_na(extract_characteristic(pheno, "race")),
    treatment = normalize_na(extract_characteristic(pheno, "treatment"))
  )

  counts_dt <- fread(count_file)
  if (!"gene_symbol" %in% colnames(counts_dt)) {
    stop("Merged GSE248493 count matrix does not contain gene_symbol column")
  }
  count_cols <- setdiff(colnames(counts_dt), "gene_symbol")
  matched_cols <- intersect(sample_map$count_sample_id, count_cols)
  extra_count_columns <- setdiff(count_cols, sample_map$count_sample_id)
  missing_count_columns <- setdiff(sample_map$count_sample_id, count_cols)
  if (length(missing_count_columns)) {
    stop("GSE248493 metadata samples missing count columns: ", paste(missing_count_columns, collapse = ", "))
  }
  if (length(extra_count_columns)) {
    message_ts("GSE248493 extra count columns excluded: ", paste(extra_count_columns, collapse = ", "))
  }

  counts_dt[, gene_symbol := clean_symbol(gene_symbol)]
  counts_dt <- counts_dt[!is.na(gene_symbol)]
  counts_dt <- counts_dt[, lapply(.SD, sum), by = gene_symbol, .SDcols = count_cols]
  symbol_map <- map_symbols_to_entrez(counts_dt$gene_symbol)
  counts_dt <- merge(counts_dt, symbol_map[mapped == TRUE, .(gene_symbol, entrez_id)], by = "gene_symbol")
  counts_dt <- counts_dt[, lapply(.SD, sum), by = entrez_id, .SDcols = matched_cols]
  counts <- as.matrix(counts_dt[, ..matched_cols])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- counts_dt$entrez_id
  sample_map <- sample_map[match(matched_cols, count_sample_id)]
  colnames(counts) <- sample_map$sample_id
  gene_dt <- data.table(
    entrez_id = rownames(counts),
    gene_symbol = map_entrez_to_symbols(rownames(counts))
  )

  list(
    gse_id = "GSE248493",
    eset = eset,
    pheno = pheno,
    sample_map = sample_map,
    expr_gene = counts,
    gene_dt = gene_dt,
    platform_id = Biobase::annotation(eset),
    data_type = "rnaseq_counts",
    expression_scale = "raw_counts_voom_logCPM",
    raw_features = nrow(counts_dt),
    gene_level_universe = nrow(counts),
    entrez_col = "mapped_from_gene_symbol",
    symbol_col = "gene_symbol",
    extra_count_columns = paste(extra_count_columns, collapse = ";"),
    missing_count_columns = paste(missing_count_columns, collapse = ";")
  )
}

sample_table_gse54837 <- function(prepared) {
  pheno <- prepared$pheno
  subject_type <- normalize_na(extract_characteristic(pheno, "subject type"))
  dt <- data.table(
    sample_id = rownames(pheno),
    geo_accession = as.character(pheno$geo_accession),
    dataset_id = "GSE54837",
    tissue = extract_characteristic(pheno, "tissue"),
    source_name_ch1 = as.character(pheno$source_name_ch1),
    subject_type = subject_type,
    age = as_numeric_clean(extract_characteristic(pheno, "age")),
    gender = factor(normalize_na(extract_characteristic(pheno, "gender"))),
    gold_stage = normalize_na(extract_characteristic(pheno, "gold stage")),
    hyb_date = normalize_na(extract_characteristic(pheno, "hyb date"))
  )
  dt[, `:=`(
    analysis_group = fifelse(
      subject_type == "COPD Subjects", "Case",
      fifelse(subject_type == "Smoker Controls", "Control", NA_character_)
    ),
    contrast_id = "GSE54837_COPD_blood_vs_smoker_control",
    contrast_label = "GSE54837 COPD blood vs smoker-control blood",
    disease_context = "COPD",
    dataset_role = "primary_copd_whole_blood_smoker_control"
  )]
  dt[, included_pre_model := !is.na(analysis_group)]
  dt[subject_type == "Non-smoker Controls", exclusion_reason := "non_smoker_control_excluded_from_primary"]
  dt[is.na(exclusion_reason) & is.na(analysis_group), exclusion_reason := "not_copd_or_smoker_control"]
  dt[is.na(exclusion_reason) & !is.na(analysis_group), exclusion_reason := "included_pre_model"]
  dt
}

sample_table_gse248493 <- function(prepared) {
  dt <- copy(prepared$sample_map)
  dt[, `:=`(
    dataset_id = "GSE248493",
    source_name_ch1 = "PBMC",
    analysis_group = fifelse(
      disease == "Chronic Obstructive Pulmonary Disease", "Case",
      fifelse(disease == "Healthy", "Control", NA_character_)
    ),
    contrast_id = "GSE248493_COPD_PBMC_RNAseq_vs_healthy",
    contrast_label = "GSE248493 COPD PBMC RNA-seq vs healthy PBMC",
    disease_context = "COPD",
    dataset_role = "supplementary_copd_pbmc_rnaseq"
  )]
  dt[, included_pre_model := !is.na(analysis_group)]
  dt[is.na(analysis_group), exclusion_reason := "not_copd_or_healthy"]
  dt[!is.na(analysis_group), exclusion_reason := "included_pre_model"]
  dt
}

sample_table_gse56766 <- function(prepared, contrast_id) {
  pheno <- prepared$pheno
  diagnosis <- normalize_na(extract_characteristic(pheno, "diagnosis"))
  stimulus <- normalize_na(extract_characteristic(pheno, "stimulus"))
  dt <- data.table(
    sample_id = rownames(pheno),
    geo_accession = as.character(pheno$geo_accession),
    dataset_id = "GSE56766",
    tissue = extract_characteristic(pheno, "tissue"),
    source_name_ch1 = as.character(pheno$source_name_ch1),
    diagnosis = diagnosis,
    stimulus = stimulus
  )

  if (identical(contrast_id, "GSE56766_COPD_unstimulated_whole_blood_vs_smoker_control")) {
    keep <- diagnosis %in% c("COPD", "Smoker") & stimulus == "Unstimulated"
    label <- "GSE56766 COPD unstimulated whole blood vs smoker-control whole blood"
    role <- "supplementary_copd_unstimulated_whole_blood"
  } else if (identical(contrast_id, "GSE56766_COPD_CD3CD28_stimulated_whole_blood_vs_smoker_control")) {
    keep <- diagnosis %in% c("COPD", "Smoker") & stimulus == "anti_CD3_anti_CD28"
    label <- "GSE56766 COPD anti-CD3/CD28-stimulated whole blood vs smoker-control whole blood"
    role <- "sensitivity_copd_tcell_stimulated_whole_blood"
  } else {
    stop("Unsupported GSE56766 contrast: ", contrast_id)
  }

  dt[, `:=`(
    analysis_group = fifelse(
      diagnosis == "COPD" & keep, "Case",
      fifelse(diagnosis == "Smoker" & keep, "Control", NA_character_)
    ),
    contrast_id = contrast_id,
    contrast_label = label,
    disease_context = "COPD",
    dataset_role = role
  )]
  dt[, included_pre_model := !is.na(analysis_group)]
  dt[is.na(stimulus), exclusion_reason := "stimulus_missing_excluded"]
  dt[diagnosis %in% c("Never Smoker", "Screen-fail"), exclusion_reason := "non_smoker_or_screen_fail_excluded"]
  dt[!is.na(stimulus) & diagnosis %in% c("COPD", "Smoker") & !keep, exclusion_reason := "other_stimulus_excluded"]
  dt[is.na(exclusion_reason) & is.na(analysis_group), exclusion_reason := "not_copd_or_smoker_control"]
  dt[is.na(exclusion_reason) & !is.na(analysis_group), exclusion_reason := "included_pre_model"]
  dt
}

sample_table_addneuromed <- function(prepared, contrast_id) {
  pheno <- prepared$pheno
  status <- normalize_na(extract_characteristic(pheno, "status"))
  included_case_control <- tolower(normalize_na(extract_characteristic(pheno, "included in case -control study")))
  dataset_id <- prepared$gse_id
  dt <- data.table(
    sample_id = rownames(pheno),
    geo_accession = as.character(pheno$geo_accession),
    dataset_id = dataset_id,
    tissue = as.character(pheno$source_name_ch1),
    status = status,
    included_case_control = included_case_control,
    ethnicity = factor(normalize_na(extract_characteristic(pheno, "ethnicity"))),
    age = as_numeric_clean(extract_characteristic(pheno, "age")),
    gender = factor(normalize_na(extract_characteristic(pheno, "gender")))
  )

  if (identical(contrast_id, paste0(dataset_id, "_AD_blood_vs_CTL_casecontrol"))) {
    keep <- status %in% c("AD", "CTL") & included_case_control == "yes"
    label <- paste0(dataset_id, " AD blood vs CTL blood, case-control subset")
    role <- if (identical(dataset_id, "GSE63060")) "primary_ad_blood_case_control_subset" else "supplementary_ad_blood_case_control_subset"
  } else if (identical(contrast_id, paste0(dataset_id, "_AD_blood_vs_CTL_all"))) {
    keep <- status %in% c("AD", "CTL")
    label <- paste0(dataset_id, " AD blood vs CTL blood, all AD/CTL samples")
    role <- "sensitivity_ad_blood_all_ad_ctl"
  } else {
    stop("Unsupported AddNeuroMed contrast: ", contrast_id)
  }

  dt[, `:=`(
    analysis_group = fifelse(
      status == "AD" & keep, "Case",
      fifelse(status == "CTL" & keep, "Control", NA_character_)
    ),
    contrast_id = contrast_id,
    contrast_label = label,
    disease_context = "AD",
    dataset_role = role
  )]
  dt[, included_pre_model := !is.na(analysis_group)]
  dt[grepl("MCI", status, ignore.case = TRUE), exclusion_reason := "mci_excluded"]
  dt[status %in% c("OTHER", "CTL to AD", "MCI to CTL"), exclusion_reason := "transition_or_other_status_excluded"]
  dt[
    grepl("_casecontrol$", contrast_id) &
      status %in% c("AD", "CTL") & included_case_control != "yes",
    exclusion_reason := "not_in_case_control_subset"
  ]
  dt[is.na(exclusion_reason) & is.na(analysis_group), exclusion_reason := "not_ad_or_control"]
  dt[is.na(exclusion_reason) & !is.na(analysis_group), exclusion_reason := "included_pre_model"]
  dt
}

build_model <- function(sample_dt, model_id, covariates) {
  dt <- copy(sample_dt[included_pre_model == TRUE])
  dt[, analysis_group := factor(analysis_group, levels = c("Control", "Case"))]

  complete_cols <- c("analysis_group", covariates)
  complete_rows <- if (length(complete_cols)) {
    stats::complete.cases(as.data.frame(dt[, ..complete_cols]))
  } else {
    rep(TRUE, nrow(dt))
  }
  dropped_covariate_missing <- sum(!complete_rows)
  dt <- dt[complete_rows]

  if (nrow(dt) < 10 || length(unique(dt$analysis_group)) != 2) {
    stop("Not enough samples after covariate filtering for ", unique(sample_dt$contrast_id), " | ", model_id)
  }
  if (min(table(dt$analysis_group)) < 5) {
    stop("One group has fewer than five samples for ", unique(sample_dt$contrast_id), " | ", model_id)
  }

  kept_covariates <- character()
  for (cov in covariates) {
    if (is.numeric(dt[[cov]])) {
      if (length(unique(dt[[cov]])) >= 2 && stats::sd(dt[[cov]]) > 0) {
        kept_covariates <- c(kept_covariates, cov)
      }
    } else {
      dt[, (cov) := factor(get(cov))]
      if (nlevels(droplevels(dt[[cov]])) >= 2) {
        kept_covariates <- c(kept_covariates, cov)
      }
    }
  }

  rhs <- paste(c("0 + analysis_group", kept_covariates), collapse = " + ")
  design <- model.matrix(as.formula(paste("~", rhs)), data = dt)
  colnames(design) <- sub("^analysis_group", "", colnames(design))
  rownames(design) <- dt$sample_id
  if (qr(design)$rank < ncol(design)) {
    stop("Design matrix is not full rank for ", unique(sample_dt$contrast_id), " | ", model_id)
  }
  contrast <- limma::makeContrasts(Case - Control, levels = design)

  list(
    sample_dt = dt,
    design = design,
    contrast = contrast,
    model_id = model_id,
    requested_covariates = covariates,
    kept_covariates = kept_covariates,
    dropped_covariate_missing = dropped_covariate_missing
  )
}

run_limma_contrast <- function(prepared, sample_dt, model_id, covariates) {
  model <- build_model(sample_dt, model_id, covariates)
  expr_sub <- prepared$expr_gene[, model$sample_dt$sample_id, drop = FALSE]
  low_count_genes_removed <- 0L
  if (identical(prepared$data_type, "rnaseq_counts")) {
    lib_size <- colSums(expr_sub)
    cpm <- t(t(expr_sub) / lib_size * 1e6)
    min_group_n <- min(table(model$sample_dt$analysis_group))
    keep <- rowSums(cpm >= 1) >= max(3L, floor(min_group_n / 2))
    low_count_genes_removed <- sum(!keep)
    expr_sub <- expr_sub[keep, , drop = FALSE]
    voom_obj <- limma::voom(expr_sub, model$design, plot = FALSE)
    fit_input <- voom_obj
    expr_for_geneset <- voom_obj$E
  } else {
    fit_input <- expr_sub
    expr_for_geneset <- expr_sub
  }
  fit <- limma::lmFit(fit_input, model$design)
  fit <- limma::contrasts.fit(fit, model$contrast)
  fit <- limma::eBayes(fit)

  gene_stats <- as.data.table(
    limma::topTable(fit, coef = 1, number = Inf, sort.by = "none"),
    keep.rownames = "entrez_id"
  )
  gene_stats <- merge(gene_stats, prepared$gene_dt, by = "entrez_id", all.x = TRUE, sort = FALSE)
  gene_stats[, `:=`(
    model_id = model_id,
    rank_metric = "limma moderated t statistic; positive = patient higher than control"
  )]
  setorder(gene_stats, -t)

  stats <- gene_stats$t
  names(stats) <- gene_stats$entrez_id
  stats <- sort(stats[is.finite(stats)], decreasing = TRUE)

  list(
    expr_sub = expr_for_geneset,
    design = model$design,
    contrast = model$contrast[, 1],
    fit = fit,
    gene_stats = gene_stats,
    rank_vector = stats,
    included_dt = model$sample_dt,
    requested_covariates = model$requested_covariates,
    kept_covariates = model$kept_covariates,
    dropped_covariate_missing = model$dropped_covariate_missing,
    low_count_genes_removed = low_count_genes_removed
  )
}

run_fgsea_axis <- function(stats, axis_entrez, pathway_name) {
  present <- intersect(axis_entrez, names(stats))
  pathways <- list(present)
  names(pathways) <- pathway_name
  fg <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    minSize = 5,
    maxSize = 500,
    eps = 0
  )
  fg <- as.data.table(fg)
  fg[, leading_edge_entrez := vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))]
  fg[, leading_edge_symbols := vapply(leadingEdge, function(x) paste(map_entrez_to_symbols(x), collapse = ";"), character(1))]
  fg[, leadingEdge := NULL]
  list(result = fg, present = present)
}



rows <- lapply(c("GSE54837", "GSE63060", "GSE63061"), function(id) {
  prepared <- prepare_microarray(id)
  sample <- if (id == "GSE54837") sample_table_gse54837(prepared) else sample_table_addneuromed(prepared, paste0(id, "_AD_blood_vs_CTL_casecontrol"))
  fit <- run_limma_contrast(prepared, sample, "age_gender_adjusted", c("age", "gender"))
  fit$gene_stats[, dataset_id := id]
  fit$gene_stats
})
fwrite(rbindlist(rows), file.path(out_dir, "gene_stats_retained_five.csv.gz"))

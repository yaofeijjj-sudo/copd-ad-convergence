# Rebuild full-gene (pre-MAD-filter) discovery matrices with the original preprocessing functions,
# verify against the saved WGCNA dat_expr, and map lung symbols to Entrez (official symbol, unique).
.libPaths(c('_workspace/r_libs_deconv_20261002','_workspace/r_libs',.libPaths()))
suppressPackageStartupMessages({library(Biobase);library(data.table);library(GEOquery);library(impute);library(limma);library(org.Hs.eg.db)})
OUT <- '03_candidate_analysis/Deconvolution_composition_adjustment_20261002'
MAX_MISSING_FRACTION <- 0.20

## ---- copied verbatim from WGCNA_figures/reproduction_materials/analysis/01_lung_wgcna/code/15_final_wgcna_copd_ad.R ----
numeric_clean <- function(x) {
  suppressWarnings(as.numeric(gsub("[^0-9.+-]", "", as.character(x))))
}
first_token <- function(x) {
  x <- as.character(x)
  x <- sub("\\s*///.*$", "", x)
  x <- sub("\\s*//.*$", "", x)
  trimws(x)
}
normalize_microarray_matrix <- function(expr_mat) {
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  if (max(expr_mat, na.rm = TRUE) > 100) {
    expr_mat <- log2(pmax(expr_mat, 1))
  }
  expr_mat <- limma::normalizeBetweenArrays(expr_mat, method = "quantile")
  keep <- rowMeans(is.na(expr_mat)) <= MAX_MISSING_FRACTION
  expr_mat <- expr_mat[keep, , drop = FALSE]
  if (anyNA(expr_mat)) {
    expr_mat <- impute::impute.knn(expr_mat)$data
  }
  expr_mat
}
collapse_probes_to_genes <- function(expr_mat, gene_id, gene_symbol = NULL,
                                     select_by = c("variance", "mean")) {
  select_by <- match.arg(select_by)
  gene_id <- first_token(gene_id)
  gene_symbol <- if (is.null(gene_symbol)) gene_id else first_token(gene_symbol)
  valid <- !is.na(gene_id) & nzchar(gene_id) & !gene_id %in% c("NA", "---")
  expr_mat <- expr_mat[valid, , drop = FALSE]
  gene_id <- gene_id[valid]
  gene_symbol <- gene_symbol[valid]
  metric <- if (select_by == "variance") apply(expr_mat, 1, var, na.rm = TRUE) else rowMeans(expr_mat, na.rm = TRUE)
  probe_table <- data.table(probe_id = rownames(expr_mat), gene_id = gene_id, gene_symbol = gene_symbol, selection_metric = metric)
  setorder(probe_table, gene_id, -selection_metric, probe_id)
  selected <- probe_table[!duplicated(gene_id)]
  out <- expr_mat[selected$probe_id, , drop = FALSE]
  rownames(out) <- selected$gene_id
  list(expr = out, map = selected[])
}
# remove_age_outliers_by_group: same rule as the original (per-group 1.5 x IQR); summary bookkeeping omitted
remove_age_outliers_by_group <- function(ages, groups) {
  keep <- rep(TRUE, length(ages))
  for (grp in sort(unique(groups))) {
    idx <- which(groups == grp); a <- ages[idx]
    q1 <- as.numeric(quantile(a, 0.25, na.rm = TRUE)); q3 <- as.numeric(quantile(a, 0.75, na.rm = TRUE)); iqr <- q3 - q1
    keep[idx[a < q1 - 1.5 * iqr | a > q3 + 1.5 * iqr]] <- FALSE
  }
  keep
}
## ---- end copied ----

## Lung GSE76925
# Offline: getGPL=FALSE, then attach the archived submitter platform tables (the tables getGEO would fetch).
PLT <- 'paper/methodology_audit_20260930/gene_identity_reconciliation_20260930/inputs'
gpl_col <- function(file,col,ids){t <- fread(file.path(PLT,file),sep='\t',colClasses='character',quote=''); stopifnot(!anyDuplicated(t$ID)); v <- t[[col]][match(ids,t$ID)]; stopifnot(!all(is.na(v))); v}
gl <- getGEO(filename='02_data/raw/bulk_lung_copd/GSE76925_series_matrix.txt.gz',GSEMatrix=TRUE,getGPL=FALSE)
el <- normalize_microarray_matrix(exprs(gl))
sym_l <- gpl_col('GPL10558_original_platform_table_20260930.tsv','Symbol',rownames(exprs(gl)))[match(rownames(el),rownames(exprs(gl)))]
ml <- collapse_probes_to_genes(el, sym_l, gene_symbol=sym_l, select_by='variance')
lung_used <- fread('WGCNA_figures/reproduction_materials/analysis/01_lung_wgcna/results/COPD_WGCNA/sample_table_used.csv')
L <- ml$expr[, lung_used$sample_id]

## Brain GSE33000
gb <- getGEO(filename='02_data/raw/bulk_brain_ad/GSE33000_series_matrix.txt.gz',GSEMatrix=TRUE,getGPL=FALSE)
pb <- pData(gb); pb <- pb[pb[['disease status:ch2']] != "Huntington's disease",,drop=FALSE]
age <- numeric_clean(pb[['age:ch2']]); dis <- data.table::fifelse(grepl('Alzheimer',pb[['disease status:ch2']],ignore.case=TRUE),1,0)
pb <- pb[remove_age_outliers_by_group(age,dis),,drop=FALSE]
eb <- normalize_microarray_matrix(exprs(gb)[,rownames(pb),drop=FALSE])
eg_b <- gpl_col('GPL4372_original_platform_table_20260930.tsv','EntrezGeneID',rownames(exprs(gb)))[match(rownames(eb),rownames(exprs(gb)))]
mb <- collapse_probes_to_genes(eb, eg_b, gene_symbol=eg_b, select_by='variance')
brain_used <- fread('01_Results_and_Methods/01_Figure1_WGCNA/analysis/02_brain_wgcna/results/AD_WGCNA_no_residual_after_age_outlier/sample_table_used.csv')
B <- mb$expr[, brain_used$sample_id]

## Verify against saved network matrices
e1 <- new.env(); load('WGCNA_figures/reproduction_materials/analysis/01_lung_wgcna/results/COPD_WGCNA/copd_wgcna_final_env.RData',envir=e1)
e2 <- new.env(); load('01_Results_and_Methods/01_Figure1_WGCNA/analysis/02_brain_wgcna/results/AD_WGCNA_no_residual_after_age_outlier/ad_wgcna_final_env.RData',envir=e2)
dl <- max(abs(t(L[colnames(e1$dat_expr),rownames(e1$dat_expr)]) - e1$dat_expr))
db <- max(abs(t(B[colnames(e2$dat_expr),rownames(e2$dat_expr)]) - e2$dat_expr))
cat('lung max|diff| vs saved dat_expr:',dl,' brain:',db,'\n')
stopifnot(dl<1e-8, db<1e-8)

## Lung symbol -> Entrez (official symbol, unique; no alias)
s2e <- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(L),keytype='SYMBOL',columns='ENTREZID')
s2e <- s2e[!is.na(s2e$ENTREZID),]; dup <- unique(s2e$SYMBOL[duplicated(s2e$SYMBOL)]); s2e <- s2e[!s2e$SYMBOL %in% dup,]
s2e <- s2e[!duplicated(s2e$ENTREZID) & !s2e$ENTREZID %in% s2e$ENTREZID[duplicated(s2e$ENTREZID)],]
Le <- L[s2e$SYMBOL,]; rownames(Le) <- s2e$ENTREZID
saveRDS(list(lung=Le,brain=B,lung_symbol_map=s2e,
             qc=list(lung_genes_total=nrow(L),lung_genes_entrez=nrow(Le),brain_genes=nrow(B),lung_maxdiff=dl,brain_maxdiff=db)),
        file.path(OUT,'source_data/bulk_full_matrices.rds'))
cat('lung genes',nrow(L),'-> Entrez',nrow(Le),'; brain genes',nrow(B),'\n')

# Follow-up of the 2 flagged AddNeuroMed pairs: within-dataset nearest-neighbour r as reference, plus age/sex/diagnosis.
.libPaths(c("_workspace/r_libs_blood_deconv_20261016", "_workspace/r_libs", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(Biobase); library(GEOquery); library(limma); library(org.Hs.eg.db); library(AnnotationDbi) })
ROOT <- normalizePath(getwd()); O <- file.path(ROOT, "03_candidate_analysis/AxisA_auxiliary_checks_20261019")
SYN <- file.path(ROOT, "GWAS_and_blood/SYN01_Entrez_downstream_20261001"); raw_dir <- file.path(ROOT, "02_data/raw/bulk_blood_external"); project_root <- file.path(O, "sandbox")
src <- parse(file.path(SYN, "code/refit_gse248493.R"))
for (e in src) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) eval(e, globalenv())
pa <- prepare_microarray("GSE63060"); pb <- prepare_microarray("GSE63061")
sa <- sample_table_addneuromed(pa, "GSE63060_AD_blood_vs_CTL_casecontrol"); sb <- sample_table_addneuromed(pb, "GSE63061_AD_blood_vs_CTL_casecontrol")
nn <- function(x) { x <- x - rowMeans(x); R <- cor(x); diag(R) <- NA; apply(R, 1, max, na.rm = TRUE) }
wa <- nn(pa$expr_gene); wb <- nn(pb$expr_gene)
cat("within-dataset nearest-neighbour r: GSE63060 median", round(median(wa), 3), "max", round(max(wa), 3), "| GSE63061 median", round(median(wb), 3), "max", round(max(wb), 3), "\n")
du <- fread(file.path(O, "results/addneuromed_cross_dataset_matches.tsv"))[flagged == TRUE]
cols <- intersect(c("sample_id", "analysis_group", "age", "gender"), names(sa))
for (i in seq_len(nrow(du))) { cat("pair", du$GSE63060[i], du$GSE63061[i], "r =", round(du$r_best[i], 3), "\n")
  print(sa[sample_id == du$GSE63060[i], ..cols]); print(sb[sample_id == du$GSE63061[i], ..cols]) }

# Build the 8 independent log-expression matrices used for Axis A replication, with identical samples and preprocessing.
# Loading code copied from AxisA_singscore_official_rerun_20261004 (02_run_specificity.R, 04_ad_reference_meta.R).
# Run from project root:
#   R_LIBS=_workspace/r_libs_affy:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_sigQC_coherence_20261014/code/01_build_matrices.R
source("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/code/00_functions.R")
suppressPackageStartupMessages(library(affy))
OUT <- "03_candidate_analysis/AxisA_sigQC_coherence_20261014"
CACHE <- "_workspace/sigqc_axisA_20261014"
for (d in c(file.path(OUT, "qa"), CACHE)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
SEL <- "GWAS_and_blood/analysis/05_disease_specificity_gwas/results/axisA_disease_context_selectivity_20260711/metadata"
RAW <- "02_data/raw/disease_specificity_tissue_20261003"
inputs <- c(
  axisA = "03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisA_confirmed_genes.csv",
  gse47460_14550 = "02_data/raw/bulk_lung_copd/GSE47460-GPL14550_series_matrix.txt.gz",
  gse47460_6480 = "02_data/raw/bulk_lung_copd/GSE47460_GPL6480/GSE47460-GPL6480_series_matrix.txt.gz",
  ph14550 = file.path(SEL, "GSE47460_GPL14550_pheno_raw.csv"), ph6480 = file.path(SEL, "GSE47460_GPL6480_pheno_raw.csv"),
  gpl14550 = "02_data/raw/platform_annotations/GPL14550_annot_table.csv", gpl6480 = "02_data/raw/platform_annotations/GPL6480_annot_table.csv",
  gse84422 = "02_data/raw/bulk_brain_ad/GSE84422-GPL96_series_matrix.txt.gz",
  ph84422 = "WGCNA_figures/reproduction_materials/analysis/02_brain_wgcna/results/AD_GSE33000_current_to_GSE84422_DLPFC_module_preservation_1000/gse84422_dlpfc_validation_metadata_and_qc.csv",
  gse118553 = "02_data/raw/bulk_brain_ad/GSE118553_series_matrix.txt.gz",
  gse313307_counts = file.path(RAW, "GSE313307_raw_count.csv.gz"), gse313307 = file.path(RAW, "GSE313307_series_matrix.txt.gz"),
  gse140841 = file.path(RAW, "GSE140841_series_matrix.txt.gz"), gse140841_tar = file.path(RAW, "GSE140841_RAW.tar"),
  gse48350 = file.path(RAW, "GSE48350_series_matrix.txt.gz"),
  gse36980 = file.path(RAW, "GSE36980_series_matrix.txt.gz"),
  gpl96 = "02_data/raw/platform_annotations/GPL96_annot_table.csv", gpl10558 = "02_data/raw/platform_annotations/GPL10558.annot.gz",
  gpl570 = "02_data/raw/platform_annotations/GPL570_annot_table.csv", gpl6244 = file.path(RAW, "GPL6244.annot.gz"),
  tx2gene = file.path(RAW, "ensembl116_tx2gene.tsv.gz"))
stopifnot(all(file.exists(P(inputs))))
fwrite(data.table(input = names(inputs), path = inputs, sha256 = sub(" .*", "", vapply(inputs, function(p) system2("shasum", c("-a", "256", shQuote(P(p))), stdout = TRUE), ""))),
       file.path(OUT, "qa/input_sha256.tsv"), sep = "\t")

# Same sample filter and log matrix as fit_comparison() in the replication branch (disease labels used only there, for filterByExpr).
prep <- function(expr, dat, covars, kind = c("array", "counts")) {
  kind <- match.arg(kind); dat <- as.data.table(dat); stopifnot(all(dat$sample %in% colnames(expr)))
  dat <- dat[complete.cases(dat[, c("disease", covars), with = FALSE])]
  e <- expr[, dat$sample, drop = FALSE]
  if (kind == "counts") {
    y <- DGEList(round(e)); keep <- filterByExpr(y, group = dat$disease); y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
    e <- cpm(y, log = TRUE, prior.count = 2)
  }
  list(expr = e, n_case = sum(dat$disease == 1), n_control = sum(dat$disease == 0))
}
mats <- list()

# 1-2 GSE47460 lung, COPD vs control donors (per-donor average), covariates age/sex/smoking
for (pl in c("14550", "6480")) {
  m <- read_matrix(inputs[[paste0("gse47460_", pl)]]); ph <- fread(P(inputs[[paste0("ph", pl)]]))
  ph <- ph[match(colnames(m), geo_accession)]; stopifnot(!anyNA(ph$geo_accession))
  m <- collapse_array(m, inputs[[paste0("gpl", pl)]], if (pl == "14550") "GENE" else "Gene ID")
  ph[, donor := sub("^(LT[0-9]{6}).*$", "\\1", title)]
  ids <- unique(ph$donor); mm <- sapply(ids, function(id) rowMeans(m[, ph$donor == id, drop = FALSE])); rownames(mm) <- rownames(m)
  pp <- ph[match(ids, donor)]
  base <- data.table(sample = ids, state = pp$`disease state:ch1`, age = suppressWarnings(as.numeric(pp$`age:ch1`)), sex = pp$`Sex:ch1`,
                     smoking = ifelse(pp$`smoker?:ch1` %in% c("", NA), "unknown", pp$`smoker?:ch1`))
  d <- base[state %in% c("Control", "Chronic Obstructive Lung Disease")][, disease := as.integer(state != "Control")]
  mats[[paste0("GSE47460-GPL", pl)]] <- prep(mm, d, c("age", "sex", "smoking"))
}
# 3 GSE84422 DLPFC
m <- read_matrix(inputs[["gse84422"]]); ph <- fread(P(inputs[["ph84422"]])); stopifnot(nrow(ph) == 57)
m <- collapse_array(m[, ph$sample_id], inputs[["gpl96"]])
mats[["GSE84422"]] <- prep(m, data.table(sample = ph$sample_id, disease = as.integer(ph$neuropathological_category != "Normal"), age = as.numeric(ph$age), sex = as.character(ph$sex_male)), c("age", "sex"))
# 4 GSE118553 frontal cortex (per individual)
md <- series_chars(inputs[["gse118553"]]); md <- md[tissue == "Frontal_Cortex" & disease.state %in% c("AD", "control")]
m <- collapse_array(read_matrix(inputs[["gse118553"]])[, md$gsm], inputs[["gpl10558"]])
ids <- unique(md$individual); mm <- sapply(ids, function(id) rowMeans(m[, md$individual == id, drop = FALSE])); rownames(mm) <- rownames(m)
mats[["GSE118553"]] <- prep(mm, md[match(ids, individual), .(sample = individual, disease = as.integer(disease.state == "AD"), age = suppressWarnings(as.numeric(age)), sex = gender)], c("age", "sex"))
# 5 GSE313307 prefrontal cortex RNA-seq
cnt <- fread(cmd = paste("gzip -cd", shQuote(P(inputs[["gse313307_counts"]]))))
mc <- as.matrix(cnt[, -1L, with = FALSE]); colnames(mc) <- trimws(colnames(mc)); ens <- sub("\\..*$", "", cnt[[1]])
ent <- unique_map(ens, "ENSEMBL"); ok <- !is.na(ent); mc <- mc[ok, ]; ent <- ent[ok]
v <- apply(log2(mc + 1), 1, mad); o <- order(ent, -v); o <- o[!duplicated(ent[o])]; mc <- mc[o, ]; rownames(mc) <- ent[o]
md <- series_chars(inputs[["gse313307"]]); md[, sample := trimws(sub(",.*$", "", title))]
mats[["GSE313307"]] <- prep(mc, md[, .(sample, disease = as.integer(ad == "1"), age = as.numeric(age), sex = gender, plate = sequencing.plate_number, infection = infection)],
                            c("age", "sex", "plate", "infection"), "counts")
# 6 GSE140841 BA9 RNA-seq (kallisto est_counts summed to Entrez)
md <- series_chars(inputs[["gse140841"]]); md <- md[tissue == "BA9" & diagnosis != "MCI"]
tmp <- file.path(tempdir(), "gse140841"); dir.create(tmp, showWarnings = FALSE); untar(P(inputs[["gse140841_tar"]]), exdir = tmp)
files <- list.files(tmp, pattern = "abundance.tsv.gz$", full.names = TRUE); fgsm <- sub("_.*$", "", basename(files))
tx <- NULL; ab <- sapply(md$gsm, function(g) { a <- fread(files[fgsm == g], select = c("target_id", "est_counts"))
  if (is.null(tx)) tx <<- sub("\\..*$", "", a$target_id) else stopifnot(identical(tx, sub("\\..*$", "", a$target_id))); a$est_counts })
t2g <- unique(fread(cmd = paste("gzip -cd", shQuote(P(inputs[["tx2gene"]]))), header = FALSE, col.names = c("tx", "gene")))
tgene <- t2g$gene[match(tx, t2g$tx)]; tent <- rep(NA_character_, length(tgene)); hasg <- !is.na(tgene)
tent[hasg] <- unique_map(tgene[hasg], "ENSEMBL"); ok <- !is.na(tent)
gc <- rowsum(ab[ok, , drop = FALSE], tent[ok]); colnames(gc) <- md$gsm
mats[["GSE140841"]] <- prep(gc, md[, .(sample = gsm, disease = as.integer(diagnosis != "Control"), age = as.numeric(age), sex = Sex, rin = as.numeric(rin.score))], c("age", "sex", "rin"), "counts")
# 7 GSE48350 superior frontal gyrus (RMA from CEL; controls aged >= youngest AD)
md <- series_chars(inputs[["gse48350"]]); md <- md[brain.region == "superior frontal gyrus"]
cels <- list.files(P(file.path(RAW, "GSE48350_SFG_CEL")), full.names = TRUE); stopifnot(setequal(sub("[_.].*$", "", basename(cels)), md$gsm))
eset <- justRMA(filenames = basename(cels), celfile.path = dirname(cels[1])); m <- exprs(eset); colnames(m) <- sub("[_.].*$", "", colnames(m))
m <- collapse_array(m, inputs[["gpl570"]])
md[, grp := fifelse(!is.na(braak.stage) | grepl("AD", individual), "AD", "control")]; md[, age := as.numeric(age..yrs.)]; minad <- min(md[grp == "AD"]$age)
mats[["GSE48350"]] <- prep(m, md[grp == "AD" | age >= minad, .(sample = gsm, disease = as.integer(grp == "AD"), age, sex = gender)], c("age", "sex"))
# 8 GSE36980 frontal cortex
md <- series_chars(inputs[["gse36980"]]); src <- system2("sh", c("-c", shQuote(sprintf("gzip -cd %s | grep '^!Sample_source_name_ch1'", shQuote(P(inputs[["gse36980"]]))))), stdout = TRUE)
md[, source := trimws(gsub('"', "", strsplit(src, "\t")[[1]][-1]))]; md <- md[tissue == "Frontal cortex"]
m <- collapse_array(read_matrix(inputs[["gse36980"]])[, md$gsm], inputs[["gpl6244"]])
mats[["GSE36980"]] <- prep(m, md[, .(sample = gsm, disease = as.integer(grepl("of AD brain", source)), age = as.numeric(age), sex = Sex)], c("age", "sex"))

# QA: counts must reproduce the replication branch
ref <- rbind(fread("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/results/axisA_specificity_results.tsv")[id %in% c("GSE47460-GPL14550 COPD", "GSE47460-GPL6480 COPD")][, id := sub(" COPD", "", id)],
             fread("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/results/AD_reference_meta.tsv")[!grepl("pooled", id)], fill = TRUE)[, .(id, ref_case = n_case, ref_control = n_control, ref_cov = axisA_covered)]
qa <- data.table(id = names(mats), tissue = c("lung", "lung", rep("frontal cortex", 6)),
                 n_case = vapply(mats, `[[`, 0L, "n_case"), n_control = vapply(mats, `[[`, 0L, "n_control"),
                 n_samples = vapply(mats, function(x) ncol(x$expr), 0L), n_genes = vapply(mats, function(x) nrow(x$expr), 0L),
                 axisA_measured = vapply(mats, function(x) sum(rownames(x$expr) %in% axisA), 0L),
                 any_nonfinite = vapply(mats, function(x) any(!is.finite(x$expr)), TRUE))
qa <- merge(qa, ref, by = "id", sort = FALSE)
qa[, matches_replication := n_case == ref_case & n_control == ref_control & axisA_measured == ref_cov]
print(qa); fwrite(qa, file.path(OUT, "qa/matrix_build_check.tsv"), sep = "\t")
stopifnot(all(qa$matches_replication), !any(qa$any_nonfinite))
saveRDS(lapply(mats, `[[`, "expr"), file.path(CACHE, "matrices.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "qa/sessionInfo_build.txt"))

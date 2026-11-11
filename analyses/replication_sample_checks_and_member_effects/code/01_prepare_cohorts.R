# Rebuild the analysed matrices of every Axis A tissue cohort, exactly as in the adopted replication
# (AxisA_singscore_official_rerun_20261004: 02_run_specificity.R for GSE47460, 04_ad_reference_meta.R for the AD cohorts),
# plus the two discovery cohorts (matrices and covariates of Figure 2d). Saves raw matrices + sample tables only;
# filtering and scoring happen in 02_ so that subsets (sensitivity analyses) are re-filtered as in the original fit_comparison().
# Run from project root:
#   R_LIBS=_workspace/r_libs_affy:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_member_direction_robustness_20261018/code/01_prepare_cohorts.R
source("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/code/00_functions.R")
suppressPackageStartupMessages(library(affy))
OUT <- "03_candidate_analysis/AxisA_member_direction_robustness_20261018"
for (d in c("source_data", "qa", "results", "figures")) dir.create(P(file.path(OUT, d)), showWarnings = FALSE)
RAW <- "02_data/raw/disease_specificity_tissue_20261003"
SEL <- "GWAS和外周血/analysis/05_disease_specificity_gwas/results/axisA_disease_context_selectivity_20260711/metadata"
coh <- list()
put <- function(id, organ, role, expr, dat, covars, kind, source, region, diagnosis, extra = list()) {
  stopifnot(all(dat$sample %in% colnames(expr)))
  coh[[id]] <<- c(list(id = id, organ = organ, role = role, expr = expr, dat = as.data.table(dat), covars = covars, kind = kind,
                       source = source, region = region, diagnosis = diagnosis), extra)
  cat(sprintf("%-18s %-6s %-10s genes=%5d samples=%4d case=%d ctrl=%d\n", id, organ, role, nrow(expr), nrow(dat), sum(dat$disease == 1), sum(dat$disease == 0)))
}

# ---------- discovery cohorts (Figure 2d inputs) ----------
b <- readRDS(P("03_candidate_analysis/Deconvolution_composition_adjustment_20261002/source_data/bulk_full_matrices.rds"))
e1 <- new.env(); load(P("WGCNA_全部图件/复现资料/analysis/01_lung_wgcna/results/COPD_WGCNA/copd_wgcna_final_env.RData"), envir = e1)
e2 <- new.env(); load(P("01_Results_and_Methods/01_Figure1_WGCNA/analysis/02_brain_wgcna/results/AD_WGCNA_no_residual_after_age_outlier/ad_wgcna_final_env.RData"), envir = e2)
cl <- as.data.table(e1$clinical_train); stopifnot(identical(cl$sample_id, colnames(b$lung)))
put("GSE76925", "lung", "discovery", b$lung, cl[, .(sample = sample_id, disease = copd_status, age = age_z, sex = sex_male, packyears = packyears_z)],
    c("age", "sex", "packyears"), "array", "Temple University / Brigham and Women's (Morrow 2017)", "lung parenchyma", "COPD (GOLD 3-4) vs smoker controls")
cb <- as.data.table(e2$clinical_train); stopifnot(identical(cb$sample_id, colnames(b$brain)))
put("GSE33000", "brain", "discovery", b$brain, cb[, .(sample = sample_id, disease = Disease_State, age = age_z, sex = sex_male)],
    c("age", "sex"), "array", "Harvard Brain Tissue Resource Center (HBTRC)", "DLPFC BA9", "AD vs non-demented control (two-colour log10 ratio)")

# ---------- lung validation: GSE47460-GPL14550 (02_run_specificity.R) ----------
m <- read_matrix("02_data/raw/bulk_lung_copd/GSE47460-GPL14550_series_matrix.txt.gz"); ph <- fread(P(file.path(SEL, "GSE47460_GPL14550_pheno_raw.csv")))
ph <- ph[match(colnames(m), geo_accession)]; stopifnot(!anyNA(ph$geo_accession))
m <- collapse_array(m, "02_data/raw/platform_annotations/GPL14550_annot_table.csv", "GENE")
ph[, donor := sub("^(LT[0-9]{6}).*$", "\\1", title)]
ids <- unique(ph$donor); mm <- sapply(ids, function(id) rowMeans(m[, ph$donor == id, drop = FALSE])); rownames(mm) <- rownames(m)
pp <- ph[match(ids, donor)]
base <- data.table(sample = ids, state = pp$`disease state:ch1`, age = suppressWarnings(as.numeric(pp$`age:ch1`)), sex = pp$`Sex:ch1`,
                   smoking = ifelse(pp$`smoker?:ch1` %in% c("", NA), "unknown", pp$`smoker?:ch1`),
                   fev1_pp = suppressWarnings(as.numeric(pp$`%predicted fev1 (pre-bd):ch1`)))
d <- base[state %in% c("Chronic Obstructive Lung Disease", "Control")][, disease := as.integer(state != "Control")]
put("GSE47460", "lung", "validation", mm, d, c("age", "sex", "smoking"), "array", "Lung Tissue Research Consortium (LGRC)", "lung", "COPD vs control (no chronic lung disease)")

# ---------- AD cohorts (04_ad_reference_meta.R) ----------
ph <- fread(P("WGCNA_全部图件/复现资料/analysis/02_brain_wgcna/results/AD_GSE33000_current_to_GSE84422_DLPFC_module_preservation_1000/gse84422_dlpfc_validation_metadata_and_qc.csv"))
m <- collapse_array(read_matrix("02_data/raw/bulk_brain_ad/GSE84422-GPL96_series_matrix.txt.gz")[, ph$sample_id], "02_data/raw/platform_annotations/GPL96_annot_table.csv")
put("GSE84422", "brain", "validation", m, data.table(sample = ph$sample_id, disease = as.integer(ph$neuropathological_category != "Normal"),
    age = as.numeric(ph$age), sex = as.character(ph$sex_male), category = ph$neuropathological_category, subject = ph$subject_id),
    c("age", "sex"), "array", "Mount Sinai/JJ Peters VA Brain Bank (MSBB)", "DLPFC BA46", "CERAD neuropathology: possible/probable/definite AD vs normal")

md <- series_chars("02_data/raw/bulk_brain_ad/GSE118553_series_matrix.txt.gz"); md <- md[tissue == "Frontal_Cortex" & disease.state %in% c("AD", "control")]
m <- collapse_array(read_matrix("02_data/raw/bulk_brain_ad/GSE118553_series_matrix.txt.gz")[, md$gsm], "02_data/raw/platform_annotations/GPL10558.annot.gz")
ids <- unique(md$individual); mm <- sapply(ids, function(id) rowMeans(m[, md$individual == id, drop = FALSE])); rownames(mm) <- rownames(m)
put("GSE118553", "brain", "validation", mm, md[match(ids, individual), .(sample = individual, disease = as.integer(disease.state == "AD"),
    age = suppressWarnings(as.numeric(age)), sex = gender)], c("age", "sex"), "array", "King's College London (submitter)", "frontal cortex", "AD vs control (neuropathology)")

cnt <- fread(cmd = paste("gzip -cd", shQuote(P(file.path(RAW, "GSE313307_raw_count.csv.gz")))))
mc <- as.matrix(cnt[, -1L, with = FALSE]); colnames(mc) <- trimws(colnames(mc)); ens <- sub("\\..*$", "", cnt[[1]])
ent <- unique_map(ens, "ENSEMBL"); ok <- !is.na(ent); mc <- mc[ok, ]; ent <- ent[ok]
v <- apply(log2(mc + 1), 1, mad); o <- order(ent, -v); o <- o[!duplicated(ent[o])]; mc <- mc[o, ]; rownames(mc) <- ent[o]
md <- series_chars(file.path(RAW, "GSE313307_series_matrix.txt.gz")); md[, sample := trimws(sub(",.*$", "", title))]
put("GSE313307", "brain", "validation", mc, md[, .(sample, disease = as.integer(ad == "1"), age = as.numeric(age), sex = gender,
    plate = sequencing.plate_number, infection = infection)], c("age", "sex", "plate", "infection"), "counts",
    "University of Exeter (submitter)", "prefrontal cortex", "AD vs control, stratified by terminal respiratory infection")

md <- series_chars(file.path(RAW, "GSE140841_series_matrix.txt.gz")); md <- md[tissue == "BA9" & diagnosis != "MCI"]
tmp <- file.path(tempdir(), "gse140841"); dir.create(tmp, showWarnings = FALSE); untar(P(file.path(RAW, "GSE140841_RAW.tar")), exdir = tmp)
files <- list.files(tmp, pattern = "abundance.tsv.gz$", full.names = TRUE); fgsm <- sub("_.*$", "", basename(files))
tx <- NULL; ab <- sapply(md$gsm, function(g) { a <- fread(files[fgsm == g], select = c("target_id", "est_counts"))
  if (is.null(tx)) tx <<- sub("\\..*$", "", a$target_id) else stopifnot(identical(tx, sub("\\..*$", "", a$target_id))); a$est_counts })
t2g <- unique(fread(cmd = paste("gzip -cd", shQuote(P(file.path(RAW, "ensembl116_tx2gene.tsv.gz")))), header = FALSE, col.names = c("tx", "gene")))
tgene <- t2g$gene[match(tx, t2g$tx)]; tent <- rep(NA_character_, length(tgene)); hasg <- !is.na(tgene)
tent[hasg] <- unique_map(tgene[hasg], "ENSEMBL"); ok <- !is.na(tent)
gc <- rowsum(ab[ok, , drop = FALSE], tent[ok]); colnames(gc) <- md$gsm
put("GSE140841", "brain", "validation", gc, md[, .(sample = gsm, disease = as.integer(diagnosis != "Control"), age = as.numeric(age), sex = Sex, rin = as.numeric(rin.score))],
    c("age", "sex", "rin"), "counts", "King's College London (submitter)", "BA9", "AD vs control (TREM2 variant carriers and non-carriers)")

md <- series_chars(file.path(RAW, "GSE48350_series_matrix.txt.gz")); md <- md[brain.region == "superior frontal gyrus"]
cels <- list.files(P(file.path(RAW, "GSE48350_SFG_CEL")), full.names = TRUE)
eset <- justRMA(filenames = basename(cels), celfile.path = dirname(cels[1])); m <- exprs(eset); colnames(m) <- sub("[_.].*$", "", colnames(m))
m <- collapse_array(m, "02_data/raw/platform_annotations/GPL570_annot_table.csv")
md[, grp := fifelse(!is.na(braak.stage) | grepl("AD", individual), "AD", "control")]; md[, age := as.numeric(age..yrs.)]; minad <- min(md[grp == "AD"]$age)
put("GSE48350", "brain", "validation", m, md[grp == "AD" | age >= minad, .(sample = gsm, disease = as.integer(grp == "AD"), age, sex = gender)],
    c("age", "sex"), "array", "UC Irvine ADRC and collaborating brain banks", "superior frontal gyrus", "AD vs control (controls aged >= youngest AD)")

md <- series_chars(file.path(RAW, "GSE36980_series_matrix.txt.gz"))
src <- system2("sh", c("-c", shQuote(sprintf("gzip -cd %s | grep '^!Sample_source_name_ch1'", shQuote(P(file.path(RAW, "GSE36980_series_matrix.txt.gz")))))), stdout = TRUE)
md[, source := trimws(gsub('"', "", strsplit(src, "\t")[[1]][-1]))]; md <- md[tissue == "Frontal cortex"]
m <- collapse_array(read_matrix(file.path(RAW, "GSE36980_series_matrix.txt.gz"))[, md$gsm], file.path(RAW, "GPL6244.annot.gz"))
put("GSE36980", "brain", "validation", m, md[, .(sample = gsm, disease = as.integer(grepl("of AD brain", source)), age = as.numeric(age), sex = Sex)],
    c("age", "sex"), "array", "Hisayama study (Kyushu University)", "frontal cortex grey matter", "AD vs non-AD (neuropathology)")

# ---------- kME of Axis A members in their discovery source modules ----------
lid <- fread(P("03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/identity/lung_all_original_features_identity.csv"), colClasses = "character")
lid <- lid[entrez_id %in% axisA & module == "red"]
kme_l <- cor(e1$dat_expr[, lid$original_gene_id], e1$MEs$MEred); kme_b <- cor(e2$dat_expr[, axisA[axisA %in% colnames(e2$dat_expr)]], e2$MEs$MEblue)
kme <- merge(data.table(entrez_id = lid$entrez_id, kME_lung_red = as.numeric(kme_l)),
             data.table(entrez_id = rownames(kme_b), kME_brain_blue = as.numeric(kme_b)), by = "entrez_id", all = TRUE)
me_dir <- data.table(module = c("lung red", "brain blue"),
  cor_with_disease = c(cor(e1$MEs$MEred, e1$clinical_train$copd_status), cor(e2$MEs$MEblue, e2$clinical_train$Disease_State)))
saveRDS(list(cohorts = coh, kme = kme, me_direction = me_dir), P(file.path(OUT, "source_data/cohorts.rds")))
fwrite(me_dir, P(file.path(OUT, "qa/source_module_eigengene_direction.tsv")), sep = "\t")
writeLines(capture.output(sessionInfo()), P(file.path(OUT, "qa/sessionInfo_01.txt")))
cat("kME rows:", nrow(kme), " lung non-NA:", sum(!is.na(kme$kME_lung_red)), " brain non-NA:", sum(!is.na(kme$kME_brain_blue)), "\n"); print(me_dir)

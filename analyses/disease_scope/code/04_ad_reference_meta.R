# AD reference: all eligible independent frontal-cortex datasets + DerSimonian-Laird random-effects pooling.
# Design: ../DESIGN_ADDENDUM_AD_reference.md ; deviations: ../DEVIATIONS.md (items 10-11).
# Run from project root:
#   R_LIBS=_workspace/r_libs_affy:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/code/04_ad_reference_meta.R
source("03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/code/00_functions.R")
suppressPackageStartupMessages({ library(metafor); library(affy) })
BR <- "03_candidate_analysis/AxisA_disease_specificity_redesign_20261003"
RAW <- "02_data/raw/disease_specificity_tissue_20261003"
inputs <- c(
  gse84422 = "02_data/raw/bulk_brain_ad/GSE84422-GPL96_series_matrix.txt.gz",
  ph84422 = "WGCNA_全部图件/复现资料/analysis/02_brain_wgcna/results/AD_GSE33000_current_to_GSE84422_DLPFC_module_preservation_1000/gse84422_dlpfc_validation_metadata_and_qc.csv",
  gse118553 = "02_data/raw/bulk_brain_ad/GSE118553_series_matrix.txt.gz",
  gse313307_counts = file.path(RAW, "GSE313307_raw_count.csv.gz"), gse313307 = file.path(RAW, "GSE313307_series_matrix.txt.gz"),
  gse140841 = file.path(RAW, "GSE140841_series_matrix.txt.gz"), gse140841_tar = file.path(RAW, "GSE140841_RAW.tar"),
  gse48350 = file.path(RAW, "GSE48350_series_matrix.txt.gz"), gse48350_cel_list = file.path(RAW, "GSE48350_SFG_cel_urls.txt"),
  gse36980 = file.path(RAW, "GSE36980_series_matrix.txt.gz"),
  gpl96 = "02_data/raw/platform_annotations/GPL96_annot_table.csv", gpl10558 = "02_data/raw/platform_annotations/GPL10558.annot.gz",
  gpl570 = "02_data/raw/platform_annotations/GPL570_annot_table.csv", gpl6244 = file.path(RAW, "GPL6244.annot.gz"),
  tx2gene = file.path(RAW, "ensembl116_tx2gene.tsv.gz"),
  addendum = file.path(BR, "DESIGN_ADDENDUM_AD_reference.md"))
stopifnot(all(file.exists(P(inputs))))
fwrite(data.table(input = names(inputs), path = inputs, sha256 = sub(" .*", "", vapply(inputs, function(p) system2("shasum", c("-a", "256", shQuote(P(p))), stdout = TRUE), ""))),
       file.path(BR, "qa/input_sha256_AD_reference.tsv"), sep = "\t")
res <- list(); flow <- list()
run <- function(seed, ...) { set.seed(seed); x <- fit_comparison(...); res[[x$id]] <<- x
  cat(sprintf("%-22s %4d/%-4d cov=%3d effect=%6.3f se=%.3f p=%.3g\n", x$id, x$n_case, x$n_control, x$axisA_covered, x$effect_SD, x$se, x$p)) }

# 1 GSE84422 DLPFC (as in the main analysis)
m <- read_matrix(inputs[["gse84422"]]); ph <- fread(P(inputs[["ph84422"]])); stopifnot(nrow(ph) == 57)
m <- collapse_array(m[, ph$sample_id], inputs[["gpl96"]])
run(1, "GSE84422", "AD_reference", "AD_dataset", "GSE84422 DLPFC", "DLPFC", "AD vs normal", m,
    data.table(sample = ph$sample_id, disease = as.integer(ph$neuropathological_category != "Normal"), age = as.numeric(ph$age), sex = as.character(ph$sex_male)), c("age", "sex"))

# 2 GSE118553 frontal cortex (as in DEVIATIONS item 10)
md <- series_chars(inputs[["gse118553"]]); md <- md[tissue == "Frontal_Cortex" & disease.state %in% c("AD", "control")]
m <- collapse_array(read_matrix(inputs[["gse118553"]])[, md$gsm], inputs[["gpl10558"]])
ids <- unique(md$individual); mm <- sapply(ids, function(id) rowMeans(m[, md$individual == id, drop = FALSE])); rownames(mm) <- rownames(m)
d <- md[match(ids, individual), .(sample = individual, disease = as.integer(disease.state == "AD"), age = suppressWarnings(as.numeric(age)), sex = gender)]
run(2, "GSE118553", "AD_reference", "AD_dataset", "GSE118553 frontal cortex", "frontal cortex", "AD vs control", mm, d, c("age", "sex"))

# 3 GSE313307 prefrontal cortex RNA-seq (infection-stratified collection)
cnt <- fread(cmd = paste("gzip -cd", shQuote(P(inputs[["gse313307_counts"]]))))
mc <- as.matrix(cnt[, -1L, with = FALSE]); colnames(mc) <- trimws(colnames(mc)); ens <- sub("\\..*$", "", cnt[[1]])
ent <- unique_map(ens, "ENSEMBL"); ok <- !is.na(ent); mc <- mc[ok, ]; ent <- ent[ok]
v <- apply(log2(mc + 1), 1, mad); o <- order(ent, -v); o <- o[!duplicated(ent[o])]; mc <- mc[o, ]; rownames(mc) <- ent[o]
md <- series_chars(inputs[["gse313307"]]); md[, sample := trimws(sub(",.*$", "", title))]; stopifnot(all(md$sample %in% colnames(mc)))
d <- md[, .(sample, disease = as.integer(ad == "1"), age = as.numeric(age), sex = gender, plate = sequencing.plate_number, infection = infection)]
run(3, "GSE313307", "AD_reference", "AD_dataset", "GSE313307 prefrontal cortex", "prefrontal cortex", "AD vs control", mc, d, c("age", "sex", "plate", "infection"), "counts")

# 4 GSE140841 BA9 RNA-seq (kallisto transcript est_counts summed to Entrez genes)
md <- series_chars(inputs[["gse140841"]]); md <- md[tissue == "BA9" & diagnosis != "MCI"]
tmp <- file.path(tempdir(), "gse140841"); dir.create(tmp, showWarnings = FALSE); untar(P(inputs[["gse140841_tar"]]), exdir = tmp)
files <- list.files(tmp, pattern = "abundance.tsv.gz$", full.names = TRUE); fgsm <- sub("_.*$", "", basename(files))
stopifnot(all(md$gsm %in% fgsm))
tx <- NULL; ab <- sapply(md$gsm, function(g) { a <- fread(files[fgsm == g], select = c("target_id", "est_counts"))
  if (is.null(tx)) tx <<- sub("\\..*$", "", a$target_id) else stopifnot(identical(tx, sub("\\..*$", "", a$target_id))); a$est_counts })
t2g <- unique(fread(cmd = paste("gzip -cd", shQuote(P(inputs[["tx2gene"]]))), header = FALSE, col.names = c("tx", "gene")))
stopifnot(!anyDuplicated(t2g$tx))
tgene <- t2g$gene[match(tx, t2g$tx)]                      # Ensembl 116 transcript -> gene
tent <- rep(NA_character_, length(tgene)); hasg <- !is.na(tgene)
tent[hasg] <- unique_map(tgene[hasg], "ENSEMBL")          # gene -> Entrez, unique only (as other RNA-seq)
ok <- !is.na(tent)
flow[["GSE140841"]] <- data.table(dataset = "GSE140841", rule = "kallisto est_counts summed by Entrez",
  transcripts = length(tx), transcripts_in_ensembl116 = sum(hasg), transcripts_mapped_entrez = sum(ok),
  median_fraction_counts_mapped = median(colSums(ab[ok, , drop = FALSE]) / colSums(ab)))
gc <- rowsum(ab[ok, , drop = FALSE], tent[ok]); colnames(gc) <- md$gsm
d <- md[, .(sample = gsm, disease = as.integer(diagnosis != "Control"), age = as.numeric(age), sex = Sex, rin = as.numeric(rin.score))]
run(4, "GSE140841", "AD_reference", "AD_dataset", "GSE140841 BA9", "BA9", "AD vs control", gc, d, c("age", "sex", "rin"), "counts")

# 5 GSE48350 superior frontal gyrus (RMA from CEL; controls restricted to age >= youngest AD case)
md <- series_chars(inputs[["gse48350"]]); md <- md[brain.region == "superior frontal gyrus"]
cels <- list.files(P(file.path(RAW, "GSE48350_SFG_CEL")), full.names = TRUE); cgsm <- sub("[_.].*$", "", basename(cels))
stopifnot(setequal(cgsm, md$gsm))
eset <- justRMA(filenames = basename(cels), celfile.path = dirname(cels[1])); m <- exprs(eset); colnames(m) <- sub("[_.].*$", "", colnames(m))
m <- collapse_array(m, inputs[["gpl570"]])
md[, grp := fifelse(!is.na(braak.stage) | grepl("AD", individual), "AD", "control")]
md[, age := as.numeric(age..yrs.)]; minad <- min(md[grp == "AD"]$age)
d <- md[grp == "AD" | age >= minad, .(sample = gsm, disease = as.integer(grp == "AD"), age, sex = gender)]
flow[["GSE48350"]] <- data.table(dataset = "GSE48350", rule = sprintf("controls aged >= %d (youngest AD)", minad), excluded_controls = md[grp == "control" & age < minad, .N])
run(5, "GSE48350", "AD_reference", "AD_dataset", "GSE48350 superior frontal gyrus", "superior frontal gyrus", "AD vs control", m, d, c("age", "sex"))

# 6 GSE36980 frontal cortex grey matter
md <- series_chars(inputs[["gse36980"]]); src <- system2("sh", c("-c", shQuote(sprintf("gzip -cd %s | grep '^!Sample_source_name_ch1'", shQuote(P(inputs[["gse36980"]]))))), stdout = TRUE)
md[, source := trimws(gsub('"', "", strsplit(src, "\t")[[1]][-1]))]; md <- md[tissue == "Frontal cortex"]
m <- collapse_array(read_matrix(inputs[["gse36980"]])[, md$gsm], inputs[["gpl6244"]])
d <- md[, .(sample = gsm, disease = as.integer(grepl("of AD brain", source)), age = as.numeric(age), sex = Sex)]
run(6, "GSE36980", "AD_reference", "AD_dataset", "GSE36980 frontal cortex", "frontal cortex", "AD vs non-AD", m, d, c("age", "sex"))

# ---------- pooling ----------
out <- rbindlist(res)
mf <- metafor::rma(yi = effect_SD, sei = se, data = out, method = "DL")
pooled <- data.table(id = "AD pooled (DL random effects)", family = "AD_reference", role = "AD_pooled", dataset = sprintf("%d datasets", nrow(out)),
  tissue = "frontal cortex", contrast = "AD vs control", n_case = sum(out$n_case), n_control = sum(out$n_control),
  effect_SD = as.numeric(mf$b), se = mf$se, ci_low = mf$ci.lb, ci_high = mf$ci.ub, p = mf$pval, I2 = mf$I2, tau2 = mf$tau2, Q = mf$QE, Q_p = mf$QEp)
loo <- rbindlist(lapply(seq_len(nrow(out)), function(i) { f <- metafor::rma(yi = effect_SD, sei = se, data = out[-i], method = "DL")
  data.table(left_out = out$id[i], effect_SD = as.numeric(f$b), ci_low = f$ci.lb, ci_high = f$ci.ub, p = f$pval, I2 = f$I2) }))
fwrite(rbind(out, pooled, fill = TRUE), file.path(BR, "results/AD_reference_meta.tsv"), sep = "\t")
fwrite(loo, file.path(BR, "results/AD_reference_leave_one_out.tsv"), sep = "\t")
fwrite(rbindlist(flow, fill = TRUE), file.path(BR, "qa/AD_reference_sample_rules.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(BR, "qa/sessionInfo_AD_reference.txt"))
print(pooled[, .(effect_SD = round(effect_SD, 3), ci_low = round(ci_low, 3), ci_high = round(ci_high, 3), p = signif(p, 3), I2 = round(I2, 1), tau2 = signif(tau2, 3), Q_p = signif(Q_p, 3))])
print(loo)

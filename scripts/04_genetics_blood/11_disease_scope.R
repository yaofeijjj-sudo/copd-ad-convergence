source("R/tissue_scores.R")
out <- "results/tissue_scope"
meta <- "results/tissue_scope/metadata"
raw <- "data/raw/tissue_scope"
annot <- "data/raw/platform_annotations"
res <- list()
add <- function(x, family, dataset, contrast) {
  res[[x$id]] <<- cbind(x[, .(id)], family = family, dataset = dataset, contrast = contrast, x[, -"id"])
}
md <- function(g) fread(file.path(meta, paste0(g, "_sample_metadata.tsv")), colClasses = "character")

# GSE47460: replicate arrays are averaged within donor.
m <- read_series_matrix("data/raw/bulk_lung_copd/GSE47460-GPL14550_series_matrix.txt.gz")
ph <- fread(file.path(meta, "GSE47460_GPL14550_pheno.csv"))
ph <- ph[match(colnames(m), geo_accession)]
stopifnot(!anyNA(ph$geo_accession))
m <- collapse_array(m, file.path(annot, "GPL14550_annot_table.csv"), "GENE")
ph[, donor := sub("^(LT[0-9]{6}).*$", "\\1", title)]
donors <- unique(ph$donor)
mm <- sapply(donors, function(d) rowMeans(m[, ph$donor == d, drop = FALSE]))
rownames(mm) <- rownames(m)
pp <- ph[match(donors, donor)]
d <- data.table(
  sample = donors, state = pp$`disease state:ch1`,
  age = suppressWarnings(as.numeric(pp$`age:ch1`)), sex = pp$`Sex:ch1`,
  smoking = ifelse(pp$`smoker?:ch1` %in% c("", NA), "unknown", pp$`smoker?:ch1`)
)
d <- d[state %in% c("Chronic Obstructive Lung Disease", "Control")][, disease := as.integer(state != "Control")]
add(fit_score("GSE47460 COPD", mm, d, c("age", "sex", "smoking")), "target", "GSE47460-GPL14550", "COPD vs control")

# GSE150910 lung RNA-seq: one shared control group for CHP, IPF and the smoking contrast.
cnt <- fread(cmd = paste("gzip -cd", file.path(raw, "GSE150910_gene-level_count_file.csv.gz")))
mc <- as.matrix(cnt[, -1L, with = FALSE])
rownames(mc) <- cnt$symbol
mc <- mc[!is.na(rownames(mc)) & nzchar(rownames(mc)), ]
ent <- unique_map(rownames(mc), "SYMBOL")
mc <- mc[!is.na(ent), ]
ent <- ent[!is.na(ent)]
if (anyDuplicated(ent)) {
  v <- apply(log2(mc + 1), 1, mad)
  o <- order(ent, -v)
  o <- o[!duplicated(ent[o])]
  mc <- mc[o, ]
  ent <- ent[o]
}
rownames(mc) <- ent
s <- md("GSE150910")
stopifnot(all(s$title %in% colnames(mc)))
base <- data.table(
  sample = s$title, grp = s$diagnosis, age = as.numeric(s$age), sex = ifelse(s$Sex == "unknown", NA, s$Sex),
  smoking = ifelse(is.na(s$ever_smoked) | s$ever_smoked == "NA", "unknown", s$ever_smoked),
  institution = s$institution, plate = s$plate
)
vs_control <- function(case) base[grp %in% c(case, "control")][, disease := as.integer(grp == case)]
cv <- c("age", "sex", "smoking", "institution", "plate")
add(fit_score("GSE150910 CHP", mc, vs_control("chp"), cv, "counts"), "lung", "GSE150910", "CHP vs control")
add(fit_score("GSE150910 IPF", mc, vs_control("ipf"), cv, "counts"), "lung", "GSE150910", "IPF vs control")
smk <- base[grp == "control" & smoking %in% c("yes", "no")][, disease := as.integer(smoking == "yes")]
add(
  fit_score("GSE150910 smoking", mc, smk, c("age", "sex", "institution", "plate"), "counts"),
  "lung", "GSE150910 controls", "ever vs never smoked (no lung disease)"
)

# Submitter DESeq2-normalized counts; the most variable Ensembl feature represents each Entrez ID.
read_norm <- function(file) {
  x <- fread(cmd = paste("gzip -cd", file.path(raw, file)))
  idcol <- if ("EnsemblID" %in% names(x)) "EnsemblID" else names(x)[1]
  keep <- setdiff(names(x), c(idcol, "symbol", "V1"))
  m <- as.matrix(x[, ..keep])
  ent <- unique_map(sub("\\..*$", "", x[[idcol]]), "ENSEMBL")
  2^keep_max_mad(log2(m[!is.na(ent), ] + 1), ent[!is.na(ent)]) - 1
}
hd <- read_norm("GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz")
s <- md("GSE64810")
stopifnot(all(s$title %in% colnames(hd)))
d <- s[, .(
  sample = title, disease = as.integer(diagnosis == "Huntington's Disease"),
  age = as.numeric(age.of.death), pmi = as.numeric(pmi), rin = as.numeric(rin)
)]
add(fit_score("GSE64810 HD", hd, d, c("age", "pmi", "rin"), "normcounts"), "brain", "GSE64810", "HD vs control")

pd <- read_norm("GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz")
s <- md("GSE68719")[, sample := sub(" .*", "", title)]
stopifnot(all(s$sample %in% colnames(pd)))
d <- s[, .(sample,
  disease = as.integer(grepl("^P_", sample)), age = as.numeric(age.at.death),
  pmi = as.numeric(post.mortem.interval..pmi.), rin = as.numeric(rna.integrity.number..rin.)
)]
add(fit_score("GSE68719 PD", pd, d, c("age", "pmi", "rin"), "normcounts"), "brain", "GSE68719", "PD vs control")

s <- md("GSE53987")[tissue == "Pre-frontal cortex (BA46)"]
m <- collapse_array(read_series_matrix(file.path(raw, "GSE53987_series_matrix.txt.gz"))[, s$gsm], file.path(annot, "GPL570_annot_table.csv"))
for (dx in c("schizophrenia", "bipolar disorder", "major depressive disorder")) {
  lab <- c(schizophrenia = "SCZ", `bipolar disorder` = "BD", `major depressive disorder` = "MDD")[[dx]]
  d <- s[
    disease.state %in% c(dx, "control"),
    .(
      sample = gsm, disease = as.integer(disease.state == dx), age = as.numeric(age), sex = gender,
      pmi = as.numeric(pmi), ph = as.numeric(ph), rin = as.numeric(rin)
    )
  ]
  add(fit_score(paste("GSE53987", lab), m, d, c("age", "sex", "pmi", "ph", "rin")), "brain", "GSE53987 BA46", paste(lab, "vs control"))
}

r <- rbindlist(res)
r[, fdr := NA_real_]
for (f in c("lung", "brain")) r[family == f, fdr := p.adjust(p, "BH")]
fwrite(r, file.path(out, "disease_scope.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_disease_scope.txt"))
print(r[, .(id, n_case, n_control, axis_a_measured, effect_SD = round(effect_SD, 2), p = signif(p, 2), fdr = signif(fdr, 2))])

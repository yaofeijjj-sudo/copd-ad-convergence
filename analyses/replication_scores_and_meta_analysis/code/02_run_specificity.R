# Axis A disease-context scope.
# Run from project root:
#   R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/code/02_run_specificity.R
suppressPackageStartupMessages({
  library(data.table); library(limma); library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db)
})
set.seed(20261003)
ROOT <- normalizePath(".")
BR <- "03_candidate_analysis/AxisA_disease_specificity_redesign_20261003"
RAW <- "02_data/raw/disease_specificity_tissue_20261003"
SEL <- "GWAS_and_blood/analysis/05_disease_specificity_gwas/results/axisA_disease_context_selectivity_20260711/metadata"
OUT <- "03_candidate_analysis/AxisA_singscore_official_rerun_20261004"
for (d in c("results", "qa")) dir.create(file.path(OUT, d), showWarnings = FALSE)
inputs <- c(
  axisA = "03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisA_confirmed_genes.csv",
  history = "03_candidate_analysis/WGCNA_Entrez_training_only_20260930/references/gene_history.tsv.gz",
  gse47460_14550 = "02_data/raw/bulk_lung_copd/GSE47460-GPL14550_series_matrix.txt.gz",
  gse47460_6480 = "02_data/raw/bulk_lung_copd/GSE47460_GPL6480/GSE47460-GPL6480_series_matrix.txt.gz",
  ph14550 = file.path(SEL, "GSE47460_GPL14550_pheno_raw.csv"),
  ph6480 = file.path(SEL, "GSE47460_GPL6480_pheno_raw.csv"),
  gse84422 = "02_data/raw/bulk_brain_ad/GSE84422-GPL96_series_matrix.txt.gz",
  ph84422 = "WGCNA_figures/reproduction_materials/analysis/02_brain_wgcna/results/AD_GSE33000_current_to_GSE84422_DLPFC_module_preservation_1000/gse84422_dlpfc_validation_metadata_and_qc.csv",
  gse33000 = "02_data/raw/bulk_brain_ad/GSE33000_series_matrix.txt.gz",
  gse118553 = "02_data/raw/bulk_brain_ad/GSE118553_series_matrix.txt.gz",
  gpl10558 = "02_data/raw/platform_annotations/GPL10558.annot.gz",
  gpl14550 = "02_data/raw/platform_annotations/GPL14550_annot_table.csv",
  gpl6480 = "02_data/raw/platform_annotations/GPL6480_annot_table.csv",
  gpl96 = "02_data/raw/platform_annotations/GPL96_annot_table.csv",
  gpl570 = "02_data/raw/platform_annotations/GPL570_annot_table.csv",
  gpl4372 = "03_candidate_analysis/WGCNA_Entrez_training_only_20260930/references/GPL4372.annot.gz",
  gse150910_counts = file.path(RAW, "GSE150910_gene-level_count_file.csv.gz"),
  gse64810_norm = file.path(RAW, "GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz"),
  gse68719_norm = file.path(RAW, "GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz"),
  gse53987 = file.path(RAW, "GSE53987_series_matrix.txt.gz"),
  gse10072 = file.path(RAW, "GSE10072_series_matrix.txt.gz"),
  meta150910 = file.path(BR, "metadata/GSE150910_sample_metadata.tsv"),
  meta64810 = file.path(BR, "metadata/GSE64810_sample_metadata.tsv"),
  meta68719 = file.path(BR, "metadata/GSE68719_sample_metadata.tsv"),
  meta53987 = file.path(BR, "metadata/GSE53987_sample_metadata.tsv"),
  meta10072 = file.path(BR, "metadata/GSE10072_sample_metadata.tsv"),
  design = file.path(BR, "DESIGN.md"))
P <- function(x) file.path(ROOT, x)
stopifnot(all(file.exists(P(inputs))))
fwrite(data.table(input = names(inputs), path = inputs,
                  sha256 = sub(" .*", "", vapply(inputs, function(p) system2("shasum", c("-a", "256", shQuote(P(p))), stdout = TRUE), ""))),
       file.path(OUT, "qa/input_sha256.tsv"), sep = "\t")

axisA <- fread(P(inputs["axisA"]), colClasses = "character")$entrez_id
stopifnot(length(axisA) == 137)
hist <- fread(cmd = paste("gzip -cd", shQuote(P(inputs["history"]))), colClasses = "character")
setnames(hist, 1, "tax_id"); hist <- hist[tax_id == "9606" & GeneID != "-"]
repl <- setNames(hist$GeneID, hist$Discontinued_GeneID)
current_id <- function(x) { for (i in 1:5) { hit <- !is.na(x) & x %in% names(repl); if (!any(hit)) break; x[hit] <- repl[x[hit]] }; x }

# ---------- readers (array rules identical to the independent-cohort branch) ----------
read_matrix <- function(file) {
  tab <- fread(cmd = paste("gzip -cd", shQuote(P(file)), "| sed -n '/^!series_matrix_table_begin/,$p' | sed '1d;/^!series_matrix_table_end/,$d'"), check.names = FALSE)
  m <- as.matrix(tab[, -1L, with = FALSE]); storage.mode(m) <- "numeric"; rownames(m) <- as.character(tab[[1L]]); m
}
keep_max_mad <- function(m, gene) {  # one feature per Entrez: largest MAD
  ok <- !is.na(gene); m <- m[ok, , drop = FALSE]; gene <- gene[ok]
  v <- apply(m, 1, mad, na.rm = TRUE); o <- order(gene, -v); o <- o[is.finite(v[o])]
  o <- o[!duplicated(gene[o])]; out <- m[o, , drop = FALSE]; rownames(out) <- gene[o]
  out[rowSums(is.finite(out)) == ncol(out), , drop = FALSE]
}
collapse_array <- function(m, annot, idcol = "Gene ID") {
  a <- if (grepl("\\.annot\\.gz$", annot)) {
    fread(cmd = paste("gzip -cd", shQuote(P(annot)), "| sed -n '/^ID\t/,$p' | grep -v '^!platform_table_end'"), colClasses = "character", select = c("ID", idcol))
  } else fread(P(annot), colClasses = "character", select = c("ID", idcol))
  setnames(a, c("probe", "gene")); a <- a[probe %in% rownames(m) & !is.na(gene) & nzchar(gene) & !grepl("/", gene)]
  keep_max_mad(m[a$probe, , drop = FALSE], current_id(a$gene))
}
unique_map <- function(keys, keytype) {  # org.Hs.eg.db, unique one-to-one key -> Entrez only, no alias
  mp <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(keys), keytype = keytype, columns = "ENTREZID"))
  mp <- as.data.table(mp)[!is.na(ENTREZID)]; setnames(mp, 1, "key")
  mp <- mp[, if (.N == 1) .SD, by = key]
  current_id(mp$ENTREZID[match(keys, mp$key)])
}
meta <- function(k) fread(P(inputs[[k]]), colClasses = "character")

# ---------- one comparison ----------
rank_matrix <- function(m) singscore::rankGenes(m)  # official singscore (2026-10-04 rerun)
fit_comparison <- function(id, family, role, dataset, tissue, contrast, expr, dat, covars, kind = c("array", "counts", "normcounts")) {
  kind <- match.arg(kind)
  dat <- as.data.table(dat); stopifnot(all(dat$sample %in% colnames(expr)))
  dat <- dat[complete.cases(dat[, c("disease", covars), with = FALSE])]
  for (v in covars) if (is.character(dat[[v]])) dat[[v]] <- factor(dat[[v]])
  covars <- covars[vapply(covars, function(v) length(unique(dat[[v]])) > 1, TRUE)]
  e <- expr[, dat$sample, drop = FALSE]
  X <- model.matrix(as.formula(paste("~ disease", if (length(covars)) paste("+", paste(covars, collapse = " + ")) else "")), dat)
  stopifnot(qr(X)$rank == ncol(X))
  if (kind == "counts") {
    y <- DGEList(round(e)); keep <- filterByExpr(y, group = dat$disease); y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
    loge <- cpm(y, log = TRUE, prior.count = 2); fit <- eBayes(lmFit(voom(y, X), X))
  } else if (kind == "normcounts") {
    keep <- filterByExpr(DGEList(round(e)), group = dat$disease); loge <- log2(e[keep, , drop = FALSE] + 1)
    fit <- eBayes(lmFit(loge, X), trend = TRUE)
  } else { loge <- e; fit <- eBayes(lmFit(loge, X)) }
  n <- nrow(loge); R <- rank_matrix(loge)
  idx <- which(rownames(loge) %in% axisA)
  score_of <- function(ii) singscore::simpleScore(R, upSet = rownames(R)[ii])$TotalScore
  eff <- function(s) { d <- copy(dat); d$score_z <- as.numeric(scale(s))
    f <- lm(as.formula(paste("score_z ~ disease", if (length(covars)) paste("+", paste(covars, collapse = " + ")) else "")), d)
    list(est = coef(f)[["disease"]], ci = confint(f)["disease", ], p = summary(f)$coefficients["disease", 4]) }
  obs <- eff(score_of(idx))
  cam <- cameraPR(fit$t[, "disease"], index = list(AxisA = idx))
  # expression-matched random sets (annotation only)
  mu <- rowMeans(loge); bins <- cut(rank(mu, ties.method = "first"), 20, labels = FALSE)
  need <- table(bins[idx]); pool <- split(setdiff(seq_len(n), idx), bins[setdiff(seq_len(n), idx)])
  null <- vapply(seq_len(1000), function(b) {
    ii <- unlist(lapply(names(need), function(k) sample(pool[[k]], need[[k]])))
    eff(score_of(ii))$est }, 0)
  emp <- min(1, 2 * min((1 + sum(null >= obs$est)) / 1001, (1 + sum(null <= obs$est)) / 1001))
  data.table(id, family, role, dataset, tissue, contrast, n_case = sum(dat$disease == 1), n_control = sum(dat$disease == 0),
             axisA_covered = length(idx), genes_ranked = n, covariates = paste(covars, collapse = ";"),
             effect_SD = obs$est, ci_low = obs$ci[[1]], ci_high = obs$ci[[2]], p = obs$p,
             camera_direction = cam$Direction, camera_p = cam$PValue,
             random_null_mean = mean(null), random_null_sd = sd(null), random_empirical_p = emp)
}
res <- list(); flow <- list()
add <- function(x) { res[[x$id]] <<- x; cat(sprintf("%-28s %4d/%-4d cov=%3d  effect=%6.3f  p=%.3g\n", x$id, x$n_case, x$n_control, x$axisA_covered, x$effect_SD, x$p)) }

# ---------- targets (reproduce the independent-cohort branch) ----------
for (pl in c("14550", "6480")) {
  m <- read_matrix(inputs[[paste0("gse47460_", pl)]]); ph <- fread(P(inputs[[paste0("ph", pl)]]))
  ph <- ph[match(colnames(m), geo_accession)]; stopifnot(!anyNA(ph$geo_accession))
  m <- collapse_array(m, inputs[[paste0("gpl", pl)]], if (pl == "14550") "GENE" else "Gene ID")
  ph[, donor := sub("^(LT[0-9]{6}).*$", "\\1", title)]
  ids <- unique(ph$donor); mm <- sapply(ids, function(id) rowMeans(m[, ph$donor == id, drop = FALSE])); rownames(mm) <- rownames(m)
  pp <- ph[match(ids, donor)]
  base <- data.table(sample = ids, state = pp$`disease state:ch1`, ild = pp$`ild subtype:ch1`,
                     age = suppressWarnings(as.numeric(pp$`age:ch1`)), sex = pp$`Sex:ch1`,
                     smoking = ifelse(pp$`smoker?:ch1` %in% c("", NA), "unknown", pp$`smoker?:ch1`))
  base[, grp := fifelse(state == "Control", "control", fifelse(state == "Chronic Obstructive Lung Disease", "COPD",
                 fifelse(state == "Interstitial lung disease" & ild == "2-UIP/IPF", "IPF", fifelse(state == "Interstitial lung disease", "otherILD", NA_character_))))]
  flow[[paste0("GSE47460-GPL", pl)]] <- base[, .N, by = grp][, dataset := paste0("GSE47460-GPL", pl)]
  cmp <- function(case, ctrl) { d <- base[grp %in% c(case, ctrl)]; d[, disease := as.integer(grp == case)]; d }
  cv <- c("age", "sex", "smoking"); ds <- paste0("GSE47460-GPL", pl)
  add(fit_comparison(paste0(ds, " COPD"), "target", if (pl == "14550") "target_primary" else "target_secondary", ds, "lung", "COPD vs control", mm, cmp("COPD", "control"), cv))
  add(fit_comparison(paste0(ds, " IPF"), "lung_secondary", "secondary", ds, "lung", "IPF/UIP vs control", mm, cmp("IPF", "control"), cv))
  add(fit_comparison(paste0(ds, " otherILD"), "lung_secondary", "secondary", ds, "lung", "other ILD vs control", mm, cmp("otherILD", "control"), cv))
  add(fit_comparison(paste0(ds, " COPD-IPF"), "lung_secondary", "secondary", ds, "lung", "COPD vs IPF/UIP", mm, cmp("COPD", "IPF"), cv))
}
m <- read_matrix(inputs[["gse84422"]]); ph <- fread(P(inputs[["ph84422"]])); stopifnot(nrow(ph) == 57)
m <- collapse_array(m[, ph$sample_id], inputs[["gpl96"]])
# GSE84422 was the frozen AD reference; replaced by GSE118553 for the figure by author decision
add(fit_comparison("GSE84422 AD", "target", "replaced_reference", "GSE84422-GPL96", "brain DLPFC", "AD vs normal", m,
                   data.table(sample = ph$sample_id, disease = as.integer(ph$neuropathological_category != "Normal"), age = as.numeric(ph$age), sex = as.character(ph$sex_male)), c("age", "sex")))

# ---------- lung comparators: GSE150910 ----------
cnt <- fread(cmd = paste("gzip -cd", shQuote(P(inputs["gse150910_counts"]))))
mc <- as.matrix(cnt[, -1L, with = FALSE]); rownames(mc) <- cnt$symbol
mc <- mc[!is.na(rownames(mc)) & nzchar(rownames(mc)), ]
ent <- unique_map(rownames(mc), "SYMBOL"); ok <- !is.na(ent)
mc <- mc[ok, ]; ent <- ent[ok]
dup <- ent[duplicated(ent)]; if (length(dup)) { v <- apply(log2(mc + 1), 1, mad); o <- order(ent, -v); o <- o[!duplicated(ent[o])]; mc <- mc[o, ]; ent <- ent[o] }
rownames(mc) <- ent
md <- meta("meta150910"); stopifnot(all(md$title %in% colnames(mc)))
base <- data.table(sample = md$title, grp = md$diagnosis, age = as.numeric(md$age), sex = ifelse(md$Sex == "unknown", NA, md$Sex),
                   smoking = ifelse(is.na(md$ever_smoked) | md$ever_smoked == "NA", "unknown", md$ever_smoked),
                   institution = md$institution, plate = md$plate)
flow[["GSE150910"]] <- base[, .N, by = grp][, dataset := "GSE150910"]
cmp <- function(case) { d <- base[grp %in% c(case, "control")]; d[, disease := as.integer(grp == case)]; d }
add(fit_comparison("GSE150910 CHP", "lung_primary", "primary", "GSE150910", "lung", "CHP vs control", mc, cmp("chp"), c("age", "sex", "smoking", "institution", "plate"), "counts"))
add(fit_comparison("GSE150910 IPF", "lung_primary", "primary", "GSE150910", "lung", "IPF vs control", mc, cmp("ipf"), c("age", "sex", "smoking", "institution", "plate"), "counts"))
sm <- base[grp == "control" & smoking %in% c("yes", "no")][, disease := as.integer(smoking == "yes")]
add(fit_comparison("GSE150910 smoking", "lung_primary", "primary", "GSE150910 controls", "lung", "ever vs never smoked (no lung disease)", mc, sm, c("age", "sex", "institution", "plate"), "counts"))

# ---------- lung secondary: GSE10072 non-tumour ----------
m <- collapse_array(read_matrix(inputs[["gse10072"]]), inputs[["gpl96"]]); md <- meta("meta10072")[source == "Normal Lung Tissue"]
stopifnot(!anyDuplicated(md$title))
d <- md[Cigarette.Smoking.Status %in% c("Current Smoker", "Never Smoked"),
        .(sample = gsm, disease = as.integer(Cigarette.Smoking.Status == "Current Smoker"), age = as.numeric(Age.at.Diagnosis), sex = Gender)]
add(fit_comparison("GSE10072 smoking", "lung_secondary", "secondary", "GSE10072 non-tumour", "lung", "current vs never smoker", m, d, c("age", "sex")))

# ---------- brain comparators ----------
read_norm <- function(k) {
  x <- fread(cmd = paste("gzip -cd", shQuote(P(inputs[k]))))
  idcol <- if ("EnsemblID" %in% names(x)) "EnsemblID" else names(x)[1]
  keep <- setdiff(names(x), c(idcol, "symbol", "V1"))
  mm <- as.matrix(x[, ..keep]); ens <- sub("\\..*$", "", x[[idcol]])
  ent <- unique_map(ens, "ENSEMBL"); ok <- !is.na(ent)
  keep_max_mad(log2(mm[ok, ] + 1), ent[ok]) -> lg
  2^lg - 1   # back to normalised-count scale for filterByExpr; ranks unchanged
}
hd <- read_norm("gse64810_norm"); md <- meta("meta64810"); stopifnot(all(md$title %in% colnames(hd)))
add(fit_comparison("GSE64810 HD", "brain_primary", "primary", "GSE64810", "brain BA9", "HD vs control", hd,
                   md[, .(sample = title, disease = as.integer(diagnosis == "Huntington's Disease"), age = as.numeric(age.of.death), pmi = as.numeric(pmi), rin = as.numeric(rin))],
                   c("age", "pmi", "rin"), "normcounts"))
pd <- read_norm("gse68719_norm"); md <- meta("meta68719"); md[, sample := sub(" .*", "", title)]; stopifnot(all(md$sample %in% colnames(pd)))
add(fit_comparison("GSE68719 PD", "brain_primary", "primary", "GSE68719", "brain BA9", "PD vs control", pd,
                   md[, .(sample, disease = as.integer(grepl("^P_", sample)), age = as.numeric(age.at.death), pmi = as.numeric(post.mortem.interval..pmi.), rin = as.numeric(rna.integrity.number..rin.))],
                   c("age", "pmi", "rin"), "normcounts"))
md <- meta("meta53987")[tissue == "Pre-frontal cortex (BA46)"]
m <- collapse_array(read_matrix(inputs[["gse53987"]])[, md$gsm], inputs[["gpl570"]])
for (dx in c("schizophrenia", "bipolar disorder", "major depressive disorder")) {
  d <- md[disease.state %in% c(dx, "control"), .(sample = gsm, disease = as.integer(disease.state == dx), age = as.numeric(age), sex = gender,
                                                   pmi = as.numeric(pmi), ph = as.numeric(ph), rin = as.numeric(rin))]
  lab <- c(schizophrenia = "SCZ", `bipolar disorder` = "BD", `major depressive disorder` = "MDD")[[dx]]
  add(fit_comparison(paste("GSE53987", lab), "brain_primary", "primary", "GSE53987 BA46", "brain BA46", paste(lab, "vs control"), m, d, c("age", "sex", "pmi", "ph", "rin")))
}
flow[["GSE53987"]] <- md[, .N, by = .(grp = disease.state)][, dataset := "GSE53987 BA46"]

# ---------- brain secondary: GSE33000 HD (controls shared with AD discovery) ----------
l <- system2("sh", c("-c", shQuote(sprintf("gzip -cd %s | grep -E '^!Sample_(geo_accession|characteristics_ch2)'", P(inputs["gse33000"])))), stdout = TRUE)
sp <- lapply(strsplit(l, "\t"), function(x) gsub('"', "", x[-1]))
gsm <- sp[[1]]; ch <- do.call(cbind, sp[-1]); field <- function(k) { j <- which(apply(ch, 2, function(x) any(startsWith(x, k)))); sub("^[^:]*: ?", "", ch[, j[1]]) }
d <- data.table(sample = gsm, status = field("disease status"), age = as.numeric(sub(" yrs", "", field("age"))), sex = field("gender"))
d <- d[status %in% c("Huntington's disease", "non-demented")][, disease := as.integer(status == "Huntington's disease")]
m <- collapse_array(read_matrix(inputs[["gse33000"]])[, d$sample], inputs[["gpl4372"]])
add(fit_comparison("GSE33000 HD", "brain_secondary", "secondary", "GSE33000", "brain PFC", "HD vs non-demented (controls shared with AD discovery)", m, d, c("age", "sex")))

# ---------- AD reference: GSE118553 frontal cortex (appended last so earlier random draws are unchanged) ----------
l <- system2("sh", c("-c", shQuote(sprintf("gzip -cd %s | grep -E '^!Sample_(geo_accession|characteristics_ch1)'", P(inputs["gse118553"])))), stdout = TRUE)
sp <- lapply(strsplit(l, "\t"), function(x) gsub('"', "", x[-1]))
ch <- do.call(cbind, sp[-1]); field <- function(k) { j <- which(apply(ch, 2, function(x) any(startsWith(x, k)))); sub("^[^:]*: ?", "", ch[, j[1]]) }
md <- data.table(gsm = sp[[1]], individual = field("individual"), tissue = field("tissue"), status = field("disease state"),
                 age = suppressWarnings(as.numeric(field("age"))), sex = field("gender"))
md <- md[tissue == "Frontal_Cortex" & status %in% c("AD", "control")]
m <- collapse_array(read_matrix(inputs[["gse118553"]])[, md$gsm], inputs[["gpl10558"]])
ids <- unique(md$individual)  # average repeated arrays within individual
mm <- sapply(ids, function(id) rowMeans(m[, md$individual == id, drop = FALSE])); rownames(mm) <- rownames(m)
d <- md[match(ids, individual), .(sample = individual, disease = as.integer(status == "AD"), age, sex)]
stopifnot(all(md[, uniqueN(status), by = individual]$V1 == 1))
flow[["GSE118553"]] <- data.table(grp = c("arrays", "individuals"), N = c(nrow(md), length(ids)), dataset = "GSE118553 frontal cortex")
add(fit_comparison("GSE118553 AD", "target", "target_primary", "GSE118553 frontal cortex", "brain frontal cortex", "AD vs control", mm, d, c("age", "sex")))

# ---------- assemble ----------
out <- rbindlist(res)
out[, fdr_family := NA_real_]
for (f in c("lung_primary", "brain_primary")) out[family == f, fdr_family := p.adjust(p, "BH")]
out[, call := fifelse(!is.na(fdr_family) & fdr_family < 0.05 & effect_SD > 0, "detected increase",
              fifelse(!is.na(fdr_family) & fdr_family < 0.05 & effect_SD < 0, "detected decrease",
              fifelse(!is.na(fdr_family), "not detected", "descriptive")))]
fwrite(out, file.path(OUT, "results/axisA_specificity_results.tsv"), sep = "\t")
fwrite(rbindlist(flow, fill = TRUE), file.path(OUT, "qa/sample_groups.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "qa/sessionInfo.txt"))
print(out[, .(id, role, n_case, n_control, axisA_covered, effect_SD = round(effect_SD, 2), ci = sprintf("%.2f to %.2f", ci_low, ci_high),
              p = signif(p, 2), fdr = signif(fdr_family, 2), call, camera = paste(camera_direction, signif(camera_p, 2)), emp = signif(random_empirical_p, 2))])

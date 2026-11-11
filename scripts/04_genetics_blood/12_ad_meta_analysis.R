source("R/tissue_scores.R")
suppressPackageStartupMessages(library(affy))
out <- "results/tissue_scope"
raw <- "data/raw/tissue_scope"
annot <- "data/raw/platform_annotations"
res <- list()
add <- function(x, dataset) res[[x$id]] <<- cbind(x[, .(id)], dataset = dataset, x[, -"id"])

ph <- fread("results/tissue_scope/metadata/GSE84422_dlpfc.csv")
m <- read_series_matrix("data/raw/bulk_brain_ad/GSE84422-GPL96_series_matrix.txt.gz")
m <- collapse_array(m[, ph$sample_id], file.path(annot, "GPL96_annot_table.csv"))
# possible, probable and definite AD are pooled against neuropathologically normal brains.
d <- data.table(
  sample = ph$sample_id, disease = as.integer(ph$neuropathological_category != "Normal"),
  age = as.numeric(ph$age), sex = as.character(ph$sex_male)
)
add(fit_score("GSE84422", m, d, c("age", "sex")), "GSE84422 DLPFC")

# GSE118553: repeated arrays are averaged within individual.
f <- "data/raw/bulk_brain_ad/GSE118553_series_matrix.txt.gz"
s <- series_chars(f)[tissue == "Frontal_Cortex" & disease.state %in% c("AD", "control")]
stopifnot(all(s[, uniqueN(disease.state), by = individual]$V1 == 1))
m <- collapse_array(read_series_matrix(f)[, s$gsm], file.path(annot, "GPL10558.annot.gz"))
ids <- unique(s$individual)
mm <- sapply(ids, function(i) rowMeans(m[, s$individual == i, drop = FALSE]))
rownames(mm) <- rownames(m)
d <- s[match(ids, individual), .(
  sample = individual, disease = as.integer(disease.state == "AD"),
  age = suppressWarnings(as.numeric(age)), sex = gender
)]
add(fit_score("GSE118553", mm, d, c("age", "sex")), "GSE118553 frontal cortex")

# GSE313307 was collected stratified by terminal respiratory infection, so infection is adjusted for.
cnt <- fread(cmd = paste("gzip -cd", file.path(raw, "GSE313307_raw_count.csv.gz")))
mc <- as.matrix(cnt[, -1L, with = FALSE])
colnames(mc) <- trimws(colnames(mc))
ent <- unique_map(sub("\\..*$", "", cnt[[1]]), "ENSEMBL")
mc <- mc[!is.na(ent), ]
ent <- ent[!is.na(ent)]
v <- apply(log2(mc + 1), 1, mad)
o <- order(ent, -v)
o <- o[!duplicated(ent[o])]
mc <- mc[o, ]
rownames(mc) <- ent[o]
s <- series_chars(file.path(raw, "GSE313307_series_matrix.txt.gz"))[, sample := trimws(sub(",.*$", "", title))]
stopifnot(all(s$sample %in% colnames(mc)))
d <- s[, .(sample,
  disease = as.integer(ad == "1"), age = as.numeric(age), sex = gender,
  plate = sequencing.plate_number, infection = infection
)]
add(fit_score("GSE313307", mc, d, c("age", "sex", "plate", "infection"), "counts"), "GSE313307 prefrontal cortex")

# GSE140841 BA9: kallisto transcript estimated counts summed to Entrez genes through Ensembl 116.
s <- series_chars(file.path(raw, "GSE140841_series_matrix.txt.gz"))[tissue == "BA9" & diagnosis != "MCI"]
tmp <- file.path(tempdir(), "GSE140841")
dir.create(tmp, showWarnings = FALSE)
untar(file.path(raw, "GSE140841_RAW.tar"), exdir = tmp)
files <- list.files(tmp, pattern = "abundance.tsv.gz$", full.names = TRUE)
fgsm <- sub("_.*$", "", basename(files))
stopifnot(all(s$gsm %in% fgsm))
tx <- NULL
ab <- sapply(s$gsm, function(g) {
  a <- fread(files[fgsm == g], select = c("target_id", "est_counts"))
  ids <- sub("\\..*$", "", a$target_id)
  if (is.null(tx)) tx <<- ids else stopifnot(identical(tx, ids))
  a$est_counts
})
t2g <- unique(fread(cmd = "gzip -cd data/reference/ensembl116_tx2gene.tsv.gz", header = FALSE, col.names = c("tx", "gene")))
stopifnot(!anyDuplicated(t2g$tx))
gene <- t2g$gene[match(tx, t2g$tx)]
tent <- rep(NA_character_, length(gene))
tent[!is.na(gene)] <- unique_map(gene[!is.na(gene)], "ENSEMBL")
gc <- rowsum(ab[!is.na(tent), , drop = FALSE], tent[!is.na(tent)])
colnames(gc) <- s$gsm
d <- s[, .(sample = gsm, disease = as.integer(diagnosis != "Control"), age = as.numeric(age), sex = Sex, rin = as.numeric(rin.score))]
add(fit_score("GSE140841", gc, d, c("age", "sex", "rin"), "counts"), "GSE140841 BA9")

# GSE48350 superior frontal gyrus, RMA from CEL files; controls younger than the youngest AD case are excluded.
s <- series_chars(file.path(raw, "GSE48350_series_matrix.txt.gz"))[brain.region == "superior frontal gyrus"]
cels <- list.files(file.path(raw, "GSE48350_SFG_CEL"), full.names = TRUE)
stopifnot(setequal(sub("[_.].*$", "", basename(cels)), s$gsm))
m <- exprs(justRMA(filenames = basename(cels), celfile.path = dirname(cels[1])))
colnames(m) <- sub("[_.].*$", "", colnames(m))
m <- collapse_array(m, file.path(annot, "GPL570_annot_table.csv"))
s[, grp := fifelse(!is.na(braak.stage) | grepl("AD", individual), "AD", "control")]
s[, age := as.numeric(age..yrs.)]
youngest <- min(s[grp == "AD"]$age)
d <- s[grp == "AD" | age >= youngest, .(sample = gsm, disease = as.integer(grp == "AD"), age, sex = gender)]
add(fit_score("GSE48350", m, d, c("age", "sex")), "GSE48350 superior frontal gyrus")

f <- file.path(raw, "GSE36980_series_matrix.txt.gz")
s <- series_chars(f)
src <- system2("sh", c("-c", shQuote(paste("gzip -cd", shQuote(f), "| grep '^!Sample_source_name_ch1'"))), stdout = TRUE)
s[, source := trimws(gsub('"', "", strsplit(src, "\t")[[1]][-1]))]
s <- s[tissue == "Frontal cortex"]
m <- collapse_array(read_series_matrix(f)[, s$gsm], file.path(annot, "GPL6244.annot.gz"))
d <- s[, .(sample = gsm, disease = as.integer(grepl("of AD brain", source)), age = as.numeric(age), sex = Sex)]
add(fit_score("GSE36980", m, d, c("age", "sex")), "GSE36980 frontal cortex")

r <- rbindlist(res)
fit <- metafor::rma(yi = effect_SD, sei = se, data = r, method = "DL")
pooled <- data.table(
  id = "pooled", dataset = sprintf("%d datasets", nrow(r)), n_case = sum(r$n_case), n_control = sum(r$n_control),
  effect_SD = as.numeric(fit$b), se = fit$se, ci_low = fit$ci.lb, ci_high = fit$ci.ub, p = fit$pval,
  I2 = fit$I2, tau2 = fit$tau2, Q = fit$QE, Q_p = fit$QEp
)
loo <- rbindlist(lapply(seq_len(nrow(r)), function(i) {
  f <- metafor::rma(yi = effect_SD, sei = se, data = r[-i], method = "DL")
  data.table(left_out = r$id[i], effect_SD = as.numeric(f$b), ci_low = f$ci.lb, ci_high = f$ci.ub, p = f$pval, I2 = f$I2)
}))
fwrite(rbind(r, pooled, fill = TRUE), file.path(out, "ad_meta_analysis.tsv"), sep = "\t")
fwrite(loo, file.path(out, "ad_leave_one_out.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_ad_meta_analysis.txt"))
print(rbind(r, pooled, fill = TRUE)[, .(id, n_case, n_control, effect_SD = round(effect_SD, 3), se = round(se, 3), p = signif(p, 3), I2)])
print(loo)

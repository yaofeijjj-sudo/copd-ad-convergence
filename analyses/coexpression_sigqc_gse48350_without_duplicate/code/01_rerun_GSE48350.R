# Rerun of the Axis A co-expression check (Table C in S2 File) for GSE48350 without GSM1176257, which the replication
# analyses have excluded since 2026-10-19 (its CEL file is byte-identical to control GSM350078 of GSE11882; Table D b).
# The 10-14 run (AxisA_sigQC_coherence_20261014) predates that exclusion and used all 41 arrays. Definitions are those of
# its 03_summarise.R (sigQC autocor_median; 1,000 random sets matched on mean expression in 20 bins; seed 20261014); the
# sample is dropped from the cached RMA matrix, as the replication branch drops it after preprocessing. The other seven
# datasets are unchanged and copied from the 10-14 results.
# Run from project root:
#   R_LIBS=_workspace/r_libs_sigqc:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_sigQC_GSE48350_dedup_20261030/code/01_rerun_GSE48350.R
suppressPackageStartupMessages(library(data.table))
SRC <- "03_candidate_analysis/AxisA_sigQC_coherence_20261014"
OUT <- "03_candidate_analysis/AxisA_sigQC_GSE48350_dedup_20261030"
for (d in c("results", "qa")) dir.create(file.path(OUT, d), showWarnings = FALSE, recursive = TRUE)
ax <- fread("03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisA_confirmed_genes.csv", colClasses = "character")
axisA <- ax$entrez_id; stopifnot(length(axisA) == 137)
excl <- fread("03_candidate_analysis/AxisA_member_direction_robustness_20261018/results/adopted_sample_exclusions_20261019.tsv")
e0 <- readRDS("_workspace/sigqc_axisA_20261014/matrices.rds")[["GSE48350"]]
drop <- excl[cohort == "GSE48350"]$sample; stopifnot(identical(drop, "GSM1176257"), drop %in% colnames(e0))
e <- e0[, setdiff(colnames(e0), drop), drop = FALSE]; stopifnot(ncol(e) == 40)

# ---- set level (as in 03_summarise.R, section 2) ----
autocor_median <- function(e) stats::median(stats::cor(t(e), method = "spearman"), na.rm = TRUE)
set.seed(20261014)
idx <- which(rownames(e) %in% axisA); obs <- autocor_median(e[idx, , drop = FALSE])
mu <- rowMeans(e); bins <- cut(rank(mu, ties.method = "first"), 20, labels = FALSE)
need <- table(bins[idx]); rest <- setdiff(seq_len(nrow(e)), idx); pool <- split(rest, bins[rest])
null <- vapply(1:1000, function(b) autocor_median(e[unlist(lapply(names(need), function(k) sample(pool[[k]], need[[k]]))), , drop = FALSE]), 0)
new <- data.table(dataset = "GSE48350", tissue = "frontal cortex", n_samples = ncol(e), axisA_measured = length(idx),
                  autocor_median = obs, null_mean = mean(null), null_q95 = unname(quantile(null, 0.95)), null_max = max(null),
                  empirical_p_one_sided = (1 + sum(null >= obs)) / 1001, coherent = obs > quantile(null, 0.95))
old <- fread(file.path(SRC, "results/set_level_coherence.tsv"))
setlev <- rbind(old[dataset != "GSE48350", names(new), with = FALSE], new)[match(old$dataset, dataset)]
fwrite(setlev, file.path(OUT, "results/set_level_coherence.tsv"), sep = "\t")

# ---- gene level (as in 03_summarise.R, section 3) ----
r <- stats::cor(t(e[idx, , drop = FALSE]), method = "spearman"); diag(r) <- NA
g48 <- data.table(entrez_id = rownames(r), GSE48350 = apply(r, 1, stats::median, na.rm = TRUE))
wide0 <- fread(file.path(SRC, "results/gene_level_coherence.tsv"), colClasses = c(entrez_id = "character"))
dcols <- c("GSE118553", "GSE140841", "GSE313307", "GSE36980", "GSE47460-GPL14550", "GSE47460-GPL6480", "GSE48350", "GSE84422")
stopifnot(all(dcols %in% names(wide0)))
wide <- merge(wide0[, c("entrez_id", "symbol", setdiff(dcols, "GSE48350")), with = FALSE], g48, by = "entrez_id", all.x = TRUE)
setcolorder(wide, c("entrez_id", "symbol", dcols)); v <- as.matrix(wide[, ..dcols])
wide[, n_measured := rowSums(!is.na(v))]; wide[, n_positive := rowSums(v > 0, na.rm = TRUE)]
wide[, median_across_datasets := apply(v, 1, stats::median, na.rm = TRUE)]
wide[, min_across_datasets := suppressWarnings(apply(v, 1, min, na.rm = TRUE))]
wide[, class := fifelse(n_measured < 6, "insufficient coverage", fifelse(n_positive == n_measured, "consistent", "inconsistent"))]
setorder(wide, -median_across_datasets)
fwrite(wide, file.path(OUT, "results/gene_level_coherence.tsv"), sep = "\t")

# ---- QA: old (41 arrays) vs new (40 arrays) ----
cmp <- merge(wide0[, .(entrez_id, old_GSE48350 = GSE48350, old_class = class)], wide[, .(entrez_id, new_GSE48350 = GSE48350, new_class = class)], by = "entrez_id")
qa <- data.table(item = c("set-level coherence", "null mean", "null 95th percentile", "null maximum", "empirical P",
                          "gene-level r (old vs new GSE48350 medians)", "max |change| in gene-level median", "members whose class changed"),
                 old = c(unlist(old[dataset == "GSE48350", .(autocor_median, null_mean, null_q95, null_max, empirical_p_one_sided)]), NA, NA, NA),
                 new = c(unlist(new[, .(autocor_median, null_mean, null_q95, null_max, empirical_p_one_sided)]),
                         cor(cmp$old_GSE48350, cmp$new_GSE48350, use = "complete.obs"),
                         max(abs(cmp$old_GSE48350 - cmp$new_GSE48350), na.rm = TRUE), sum(cmp$old_class != cmp$new_class)))
print(qa); fwrite(qa, file.path(OUT, "qa/old_vs_new_GSE48350.tsv"), sep = "\t")
print(wide[, .N, by = class])
writeLines(capture.output(sessionInfo()), file.path(OUT, "qa/sessionInfo.txt"))

# Advice 5.1: (a) confirm every module set used in the 128 Fisher tests is restricted to the same 3,623-gene background;
# (b) repeat all 128 one-sided Fisher tests (stats::fisher.test, alternative = "greater", BH over 128) with the
# grey-inclusive background of 4,832 shared confirmed IDs. Module definitions are unchanged; only the universe changes.
# Run from project root: _workspace/envs/r/bin/Rscript 03_candidate_analysis/Methods_checks_20261019/code/01_fisher_background.R
suppressPackageStartupMessages(library(data.table))
IN <- "03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross"
OUT <- "03_candidate_analysis/Methods_checks_20261019"
g <- fread(file.path(IN, "shared_all_module_gene_ids.csv"), colClasses = "character")
orig <- fread(file.path(IN, "all128_module_pairs_fisher_BH.csv"))
stopifnot(nrow(g) == 4832, uniqueN(g$entrez_id) == 4832)
run <- function(bg, label) {
  L <- setdiff(unique(g$module_lung), "grey"); B <- setdiff(unique(g$module_brain), "grey"); N <- nrow(bg)
  r <- rbindlist(lapply(L, function(l) rbindlist(lapply(B, function(b) {
    a <- sum(bg$module_lung == l & bg$module_brain == b); ll <- sum(bg$module_lung == l); bb <- sum(bg$module_brain == b)
    ft <- fisher.test(matrix(c(a, bb - a, ll - a, N - ll - bb + a), 2, byrow = TRUE), alternative = "greater")
    data.table(lung_module = l, brain_module = b, overlap = a, lung_genes = ll, brain_genes = bb, background = N, odds_ratio = unname(ft$estimate), p_value = ft$p.value) }))))
  r[, fdr := p.adjust(p_value, "BH")][, background_label := label][]
}
ng <- run(g[module_lung != "grey" & module_brain != "grey"], "3,623 non-grey (original)")
gr <- run(g, "4,832 grey-inclusive")
# (a) reproduce the original table exactly
chk <- merge(ng, orig[, .(lung_module, brain_module, o_overlap = overlap, o_l = lung_genes, o_b = brain_genes, o_p = p_value, o_fdr = fdr)], by = c("lung_module", "brain_module"))
stopifnot(nrow(chk) == 128, all(chk$overlap == chk$o_overlap), all(chk$lung_genes == chk$o_l), all(chk$brain_genes == chk$o_b), max(abs(chk$p_value - chk$o_p) / pmax(chk$o_p, 1e-300)) < 1e-8)
both <- merge(ng, gr, by = c("lung_module", "brain_module"), suffixes = c("_nongrey", "_grey"))
both[, sig_nongrey := fdr_nongrey < 0.05][, sig_grey := fdr_grey < 0.05]
fwrite(gr, file.path(OUT, "results/fisher_128_pairs_grey_inclusive.tsv"), sep = "\t")
fwrite(both, file.path(OUT, "results/fisher_128_pairs_background_comparison.tsv"), sep = "\t")
cat("significant pairs (BH<0.05): non-grey", sum(both$sig_nongrey), " grey-inclusive", sum(both$sig_grey), " both", sum(both$sig_nongrey & both$sig_grey), "\n")
print(both[(lung_module == "red" & brain_module == "blue") | (lung_module == "turquoise" & brain_module == "turquoise"),
  .(lung_module, brain_module, overlap_nongrey, lung_genes_nongrey, brain_genes_nongrey, OR_ng = round(odds_ratio_nongrey, 2), fdr_ng = signif(fdr_nongrey, 3),
    overlap_grey, lung_genes_grey, brain_genes_grey, OR_g = round(odds_ratio_grey, 2), fdr_g = signif(fdr_grey, 3))])
print(both[sig_nongrey != sig_grey, .(lung_module, brain_module, fdr_nongrey = signif(fdr_nongrey, 3), fdr_grey = signif(fdr_grey, 3))])
cat("rank of red x blue by P: non-grey", rank(both$p_value_nongrey)[both$lung_module == "red" & both$brain_module == "blue"], " grey", rank(both$p_value_grey)[both$lung_module == "red" & both$brain_module == "blue"], "\n")

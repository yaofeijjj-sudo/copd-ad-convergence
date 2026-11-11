source("R/runtime.R")
source("scripts/05_single_cell/common.R")
LUNG <- file.path(ROOT, "results/single_cell/results/lung")
OUT <- file.path(ROOT, "results/single_cell/results/models")
x <- fread(file.path(LUNG, "objects/axisA_cell_scores.tsv.gz"), select = c("gsm", "cell_id", "subject_id", "broad_celltype", "AxisA_UCell", "celltypist_label"))
stopifnot(nrow(x) == 211464, uniqueN(x$cell_id) == nrow(x), all(is.finite(x$AxisA_UCell)))
x <- x[broad_celltype == "Myeloid"]
stopifnot(nrow(x) == 116558)
map <- data.table(fine_label = c("Alveolar macrophages", "Alveolar Mph proliferating", "Alveolar Mph CCL3+", "Alveolar Mph MT-positive", "Monocyte-derived Mph", "Interstitial Mph perivascular", "Classical monocytes", "Non-classical monocytes", "DC1", "DC2", "Migratory DCs", "Plasmacytoid DCs", "Mast cells"), group = c(rep("Alveolar macrophages", 4), "Monocyte-derived macrophages", "Perivascular interstitial macrophages", "Classical monocytes", "Non-classical monocytes", rep("Conventional DCs", 3), "Plasmacytoid DCs", "Mast cells"))
stopifnot(setequal(unique(x$celltypist_label), map$fine_label))
x[, annotation := map$group[match(celltypist_label, map$fine_label)]]
write_tsv(map, "lung_myeloid_label_map.tsv")
lib <- x[, .(n_cells = .N, mean_ucell = mean(AxisA_UCell)), by = .(gsm, donor = subject_id, annotation)]
g <- lib[, .(n_cells = sum(n_cells), n_libraries = .N, mean_ucell = mean(mean_ucell)), by = .(donor, annotation)]
write_tsv(lib, "lung_myeloid_library_profiles.tsv")
write_tsv(g, "lung_myeloid_donor_profiles.tsv")
cov <- g[n_cells >= 20, .(n_donors = uniqueN(donor), n_cells = sum(n_cells)), by = annotation]
types <- sort(cov[n_donors >= 5]$annotation)
stopifnot(length(types) == 8)
model <- fit_model(g, types)
res <- rbindlist(lapply(types, function(target) {
  w <- setNames(rep(-1 / (length(types) - 1), length(types)), types)
  w[target] <- 1
  contrast(model, w, "GSE302339", target, "myeloid_subgroup")$row
}))
res[, BH_FDR := p.adjust(p, "BH")]
res <- merge(res, cov, by.x = "contrast", by.y = "annotation", suffixes = c("_model", "_target"))
write_tsv(res, "lung_myeloid_CR2.tsv")
write_tsv(cov, "lung_myeloid_group_coverage.tsv")
fl <- x[, .(n_cells = .N, mean_ucell = mean(AxisA_UCell)), by = .(gsm, donor = subject_id, annotation = celltypist_label)]
fd <- fl[, .(n_cells = sum(n_cells), n_libraries = .N, mean_ucell = mean(mean_ucell)), by = .(donor, annotation)]
write_tsv(fd, "lung_myeloid_fine_donor_profiles.tsv")
fl[, group := map$group[match(annotation, map$fine_label)]]
v <- fl[, .(mean_ucell = weighted.mean(mean_ucell, n_cells), n_cells = sum(n_cells)), by = .(gsm, donor, annotation = group)]
v <- v[, .(check_mean = mean(mean_ucell), check_n = sum(n_cells)), by = .(donor, annotation)]
j <- merge(g, v, by = c("donor", "annotation"))
stopifnot(nrow(j) == nrow(g), all(j$n_cells == j$check_n), max(abs(j$mean_ucell - j$check_mean)) < 1e-12)
qa <- data.table(check = c("all myeloid cells retained", "independent library-aware grouping", "model full rank", "all eight prespecified groups eligible"), value = c(sum(g$n_cells), max(abs(j$mean_ucell - j$check_mean)), model$fit$rank, length(types)), pass = TRUE)
write_tsv(qa, "lung_myeloid_numeric_QA.tsv")
saveRDS(model, file.path(OUT, "objects/lung_myeloid_model.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "qa/sessionInfo_lung.txt"))
if ("lung-only" %in% commandArgs(trailingOnly = TRUE)) quit(save = "no", status = 0)
b <- fread(file.path(OUT, "tables/brain_subtype_CR2.tsv"))
bm <- b[, .(original_label = contrast, display_label = gsub("_", " ", gsub("-SEAAD", "", contrast, fixed = TRUE), fixed = TRUE))]
bm[, display_label := gsub("Micro-PVM 2 3", "Micro-PVM 2.3", display_label, fixed = TRUE)]
bm[, display_label := gsub("Micro-PVM 2 1", "Micro-PVM 2.1", display_label, fixed = TRUE)]
bm[original_label == "Lymphocyte", display_label := "Lymphocytes"]
stopifnot(uniqueN(bm$display_label) == nrow(bm))
write_tsv(bm, "brain_subtype_display_map.tsv")
print(res[, .(contrast, effect, ci_low, ci_high, BH_FDR, n_donors_target)])
print(qa)

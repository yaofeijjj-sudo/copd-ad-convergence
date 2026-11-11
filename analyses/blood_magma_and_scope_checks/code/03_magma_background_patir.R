# Advice 3.2 (DESIGN.md): (a) four MAGMA tests reproduced, then rerun with the background restricted to the 4,832 shared
# network-eligible genes (MAGMA v1.10 --settings gene-include); (b) AD test conditional on Patir 2019 core microglia (249 genes).
# Run from project root:
#   LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs Rscript 03_candidate_analysis/AxisA_auxiliary_checks_20261019/code/03_magma_background_patir.R
suppressPackageStartupMessages({ library(data.table); library(AnnotationDbi); library(org.Hs.eg.db) })
O <- "03_candidate_analysis/AxisA_auxiliary_checks_20261019"; for (d in c("magma", "inputs", "logs")) dir.create(file.path(O, d), FALSE)
magma <- "_workspace/project_tools/external_tools/magma_v1.10_official/magma_v1.10_mac/magma"
raw <- c(AD = "GWAS和外周血/SYN01_Entrez_downstream_20261001/results/magma/AD_total_N.genes.raw",
         COPD = "GWAS和外周血/analysis/05_disease_specificity_gwas/results/current_axes_minimal_ad_copd_magma_20260802/magma/GCST90244098.genic_only.genes.raw")
sets_f <- "GWAS和外周血/SYN01_Entrez_downstream_20261001/results/magma_inputs/locked_axes.sets"
frozen <- fread("GWAS和外周血/SYN01_Entrez_downstream_20261001/results/magma_four_tests.csv")
bg <- fread("03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/shared_all_module_gene_ids.csv", colClasses = "character")$entrez_id
stopifnot(length(bg) == 4832); bg_f <- file.path(O, "inputs/background_4832_shared_ids.txt"); writeLines(bg, bg_f)
sets <- lapply(strsplit(readLines(sets_f), "\t"), function(x) x[-1]); names(sets) <- sapply(strsplit(readLines(sets_f), "\t"), `[`, 1)
stopifnot(all(sets$AxisA %in% bg), all(sets$AxisB %in% bg))
run <- function(args, log) { out <- system2(magma, args, stdout = TRUE, stderr = TRUE); writeLines(c(paste(args, collapse = " "), out), file.path(O, "logs", log)); stopifnot(is.null(attr(out, "status"))) }
gsa <- function(f) fread(cmd = paste("grep -v '^#'", shQuote(f)))

res <- list()
for (dz in names(raw)) for (b in c("all genes (original)", "4,832 network-eligible")) {
  tag <- paste0(dz, if (b == "all genes (original)") "_orig" else "_bg4832")
  args <- c("--gene-results", raw[[dz]], "--set-annot", sets_f, "--out", file.path(O, "magma", tag))
  if (b != "all genes (original)") args <- c(args, "--settings", paste0("gene-include=", bg_f))
  run(args, paste0(tag, ".log")); g <- gsa(file.path(O, "magma", paste0(tag, ".gsa.out")))
  g[, `:=`(disease = dz, background = b)]; res[[tag]] <- g }
res <- rbindlist(res, fill = TRUE)
res[, holm_4tests := p.adjust(P, "holm"), by = background]
res[, `:=`(ci_low = BETA - 1.96 * SE, ci_high = BETA + 1.96 * SE)]
chk <- merge(res[background == "all genes (original)", .(VARIABLE, disease, BETA, P)], frozen[, .(VARIABLE, disease, BETA_f = BETA, P_f = P)], by = c("VARIABLE", "disease"))
print(chk); stopifnot(nrow(chk) == 4, max(abs(chk$BETA - chk$BETA_f)) < 1e-4, max(abs(chk$P - chk$P_f) / chk$P_f) < 1e-3)
fwrite(res, file.path(O, "results/magma_background_sensitivity.tsv"), sep = "\t")
print(res[, .(disease, VARIABLE, background, NGENES, BETA = round(BETA, 3), ci = sprintf("%.3f to %.3f", ci_low, ci_high), P = signif(P, 3), holm_4tests = signif(holm_4tests, 3))])

# ---- Patir 2019 core microglia (Table S3, "Core" = 249 genes), AD only, same model as the 10-16 BRETIGEA analysis ----
core <- readLines(file.path(O, "source_data/patir2019_core_symbols.txt")); stopifnot(length(core) == 249)   # extracted from Table S3 column "Core signature" == "Core" (code/03b_extract_patir.py)
m <- AnnotationDbi::select(org.Hs.eg.db, keys = core, keytype = "SYMBOL", columns = "ENTREZID"); m <- m[!is.na(m$ENTREZID), ]
m <- m[!(m$SYMBOL %in% m$SYMBOL[duplicated(m$SYMBOL)]), ]; patir <- unique(m$ENTREZID)
analysed <- sub("^\\s*(\\S+).*$", "\\1", grep("^#", readLines(raw[["AD"]]), invert = TRUE, value = TRUE))
man <- data.table(n_core_symbols = 249, n_entrez = length(patir), n_in_gene_analysis = sum(patir %in% analysed),
  overlap_AxisA = sum(sets$AxisA %in% patir), overlap_AxisA_analysed = sum(sets$AxisA %in% patir & sets$AxisA %in% analysed),
  unmapped_symbols = paste(setdiff(core, m$SYMBOL), collapse = ";"))
fwrite(man, file.path(O, "results/patir_core_manifest.tsv"), sep = "\t"); print(man)
sf <- file.path(O, "inputs/axes_plus_patir.sets")
writeLines(c(paste(c("AxisA", sets$AxisA), collapse = "\t"), paste(c("AxisB", sets$AxisB), collapse = "\t"),
             paste(c("PatirCore", patir), collapse = "\t"), paste(c("AxisA_nonPatir", setdiff(sets$AxisA, patir)), collapse = "\t")), sf)
run(c("--gene-results", raw[["AD"]], "--set-annot", sf, "--out", file.path(O, "magma", "AD_patir_uncond")), "AD_patir_uncond.log")
run(c("--gene-results", raw[["AD"]], "--set-annot", sf, "--model", "condition=PatirCore", "analyse=list,AxisA,AxisB", "--out", file.path(O, "magma", "AD_patir_cond")), "AD_patir_cond.log")
pt <- rbind(cbind(gsa(file.path(O, "magma", "AD_patir_uncond.gsa.out")), model = "unconditional"),
            cbind(gsa(file.path(O, "magma", "AD_patir_cond.gsa.out")), model = "conditional on PatirCore"), fill = TRUE)
pt[VARIABLE %in% c("AxisA", "AxisB"), holm_two_axes := p.adjust(P, "holm"), by = model]
fwrite(pt, file.path(O, "results/magma_patir_conditional.tsv"), sep = "\t")
print(pt[, .(model, VARIABLE, NGENES, BETA = round(BETA, 3), SE = round(SE, 3), P = signif(P, 3), holm_two_axes = signif(holm_two_axes, 3))])

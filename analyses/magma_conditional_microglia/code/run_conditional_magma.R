# R7 conditional MAGMA. Run from the project root:
#   LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs_deconv_20261002:_workspace/r_libs Rscript 03_candidate_analysis/AxisA_MAGMA_conditional_microglia_20261016/code/run_conditional_magma.R
suppressPackageStartupMessages({library(BRETIGEA); library(AnnotationDbi); library(org.Hs.eg.db); library(data.table)})
O <- "03_candidate_analysis/AxisA_MAGMA_conditional_microglia_20261016"
for (d in c("inputs", "magma", "results", "logs")) dir.create(file.path(O, d), showWarnings = FALSE)
magma <- "_workspace/project_tools/external_tools/magma_v1.10_official/magma_v1.10_mac/magma"
genes_raw <- "GWAS_and_blood/SYN01_Entrez_downstream_20261001/results/magma/AD_total_N.genes.raw"
base_sets <- "GWAS_and_blood/SYN01_Entrez_downstream_20261001/results/magma_inputs/locked_axes.sets"
stopifnot(as.character(packageVersion("BRETIGEA")) == "1.0.3", as.character(packageVersion("org.Hs.eg.db")) == "3.22.0")

sets <- lapply(strsplit(readLines(base_sets), "\t"), function(x) x[-1]); names(sets) <- sapply(strsplit(readLines(base_sets), "\t"), `[`, 1)
stopifnot(length(sets$AxisA) == 137, length(sets$AxisB) == 434)
raw_lines <- readLines(genes_raw); raw_lines <- raw_lines[!startsWith(raw_lines, "#")]
analysed <- sub("^\\s*(\\S+).*$", "\\1", raw_lines)   # first whitespace-delimited field = Entrez ID (rows have variable length)
stopifnot(length(analysed) == 18224)

data(markers_df_human_brain)
mic <- markers_df_human_brain$markers[markers_df_human_brain$cell == "mic"]
stopifnot(length(mic) == 1000)
map_sym <- function(sym) {
  m <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(sym), keytype = "SYMBOL", columns = "ENTREZID")
  m <- m[!is.na(m$ENTREZID), ]; m <- m[!(m$SYMBOL %in% m$SYMBOL[duplicated(m$SYMBOL)]), ]
  ids <- m$ENTREZID[match(sym, m$SYMBOL)]; ids[!is.na(ids)]
}
ref <- list(top250 = map_sym(mic[1:250]), all1000 = map_sym(mic))
manifest <- rbindlist(lapply(names(ref), function(k) data.table(reference = k, n_symbols = if (k == "top250") 250 else 1000,
  n_entrez = length(ref[[k]]), n_in_gene_analysis = sum(ref[[k]] %in% analysed),
  overlap_AxisA = sum(sets$AxisA %in% ref[[k]]), overlap_AxisA_analysed = sum(sets$AxisA %in% ref[[k]] & sets$AxisA %in% analysed),
  overlap_AxisB = sum(sets$AxisB %in% ref[[k]]))))
fwrite(manifest, file.path(O, "results/reference_set_manifest.csv"))

run <- function(args, log) {
  out <- system2(magma, args, stdout = TRUE, stderr = TRUE)
  writeLines(c(paste(args, collapse = " "), out), file.path(O, "logs", log))
  stopifnot(is.null(attr(out, "status")))
}
read_gsa <- function(f) fread(cmd = paste("grep -v '^#'", shQuote(f)))
res <- list()
for (k in names(ref)) {
  sf <- file.path(O, "inputs", paste0("axes_plus_microglia_", k, ".sets"))
  writeLines(c(paste(c("AxisA", sets$AxisA), collapse = "\t"), paste(c("AxisB", sets$AxisB), collapse = "\t"),
               paste(c("MicrogliaRef", ref[[k]]), collapse = "\t"),
               paste(c("AxisA_nonRef", setdiff(sets$AxisA, ref[[k]])), collapse = "\t")), sf)
  # unconditional (all four sets, each alone)
  run(c("--gene-results", genes_raw, "--set-annot", sf, "--out", file.path(O, "magma", paste0("AD_uncond_", k))), paste0("uncond_", k, ".log"))
  u <- read_gsa(file.path(O, "magma", paste0("AD_uncond_", k, ".gsa.out")))
  u[, `:=`(reference = k, model = "unconditional")]
  # conditional on MicrogliaRef, analysing the two axes
  run(c("--gene-results", genes_raw, "--set-annot", sf, "--model", "condition=MicrogliaRef", "analyse=list,AxisA,AxisB",
        "--out", file.path(O, "magma", paste0("AD_cond_", k))), paste0("cond_", k, ".log"))
  cnd <- read_gsa(file.path(O, "magma", paste0("AD_cond_", k, ".gsa.out")))
  cnd[, `:=`(reference = k, model = "conditional on MicrogliaRef")]
  res[[k]] <- rbind(u, cnd, fill = TRUE)
}
res <- rbindlist(res, fill = TRUE)
res[model == "conditional on MicrogliaRef" & VARIABLE %in% c("AxisA", "AxisB"), holm_P_two_axes := p.adjust(P, "holm"), by = reference]
res[model == "unconditional" & VARIABLE %in% c("AxisA", "AxisB"), holm_P_two_axes := p.adjust(P, "holm"), by = reference]
fwrite(res, file.path(O, "results/conditional_magma_results.csv"))
writeLines(capture.output(sessionInfo()), file.path(O, "logs/sessionInfo.txt"))
print(manifest); print(res[, .(reference, model, VARIABLE, NGENES, BETA, SE, P, holm_P_two_axes)])

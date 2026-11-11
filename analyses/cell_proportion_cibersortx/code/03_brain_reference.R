# SEA-AD DFC reference -> CIBERSORTx text files (version A all genes; version B without Axis A), same gene mapping as the lung reference.
# The website limit is 50 million data points; if the A version exceeds it, nuclei per type are reduced to 200 (pre-specified fallback).
# Run from project root:
#   R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018/code/03_brain_reference.R
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(AnnotationDbi); library(org.Hs.eg.db) })
OUT <- "03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018"
X <- as(readMM(file.path(OUT, "source_data/brain_ref_counts.mtx")), "CsparseMatrix")
genes <- fread(file.path(OUT, "source_data/brain_ref_genes.tsv"))$ensembl; cells <- fread(file.path(OUT, "source_data/brain_ref_cells.tsv"))
stopifnot(nrow(X) == length(genes), ncol(X) == nrow(cells))
mp <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", columns = "ENTREZID"))[!is.na(ENTREZID)]
mp <- mp[!ENSEMBL %in% ENSEMBL[duplicated(ENSEMBL)] & !ENTREZID %in% ENTREZID[duplicated(ENTREZID)]]
mp[, SYMBOL := mapIds(org.Hs.eg.db, ENTREZID, "SYMBOL", "ENTREZID")]; mp <- mp[!is.na(SYMBOL) & !SYMBOL %in% SYMBOL[duplicated(SYMBOL)]]
X <- X[match(mp$ENSEMBL, genes), ]; rownames(X) <- mp$SYMBOL; ent <- mp$ENTREZID
pick <- function(n) { set.seed(20261018); cells[, k := seq_len(.N)]; s <- cells[, .SD[sort(sample(.N, min(n, .N)))], by = type]$k; s }
cols <- seq_len(ncol(X)); keep <- rowSums(X[, cols]) > 0
if (sum(keep) * length(cols) > 5e7) { cols <- pick(200); keep <- rowSums(X[, cols]) > 0; message("fallback: 200 nuclei per type") }
X <- X[keep, cols]; ent <- ent[keep]; labs <- cells$type[cols]
axisA <- fread("03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisA_confirmed_genes.csv", colClasses = "character")$entrez_id
write_ref <- function(M, file) { con <- file(file, "w"); writeLines(paste(c("GeneSymbol", labs), collapse = "\t"), con)
  M <- as.matrix(M); storage.mode(M) <- "integer"; write.table(M, con, sep = "\t", quote = FALSE, col.names = FALSE); close(con) }
write_ref(X, file.path(OUT, "upload/brain_scref_A.txt")); write_ref(X[!ent %in% axisA, ], file.path(OUT, "upload/brain_scref_B.txt"))
saveRDS(list(genes_symbol = rownames(X), entrez = ent, labels = labs, cells = cells[cols]), file.path(OUT, "source_data/brain_reference_meta.rds"))
fwrite(as.data.table(table(type = labs)), file.path(OUT, "qa/brain_reference_cells_per_type.tsv"), sep = "\t")
cat("nuclei", ncol(X), " genes A", nrow(X), " genes B", sum(!ent %in% axisA), " Axis A removed", sum(ent %in% axisA), " datapoints A", ncol(X) * nrow(X), "\n"); print(table(labs))

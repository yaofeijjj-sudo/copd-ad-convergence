# GSE302339 single-cell reference for CIBERSORTx (DESIGN.md): 300 QC cells per broad type, balanced across donors,
# raw UMI counts, gene symbols; version A (all genes) and version B (137 Axis A genes removed).
# Run from project root:
#   R_LIBS=_workspace/r_libs_deconv_20261002:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018/code/01_lung_reference.R
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(hdf5r); library(AnnotationDbi); library(org.Hs.eg.db) })
OUT <- "03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018"
for (d in c("upload", "qa", "results", "source_data")) dir.create(file.path(OUT, d), showWarnings = FALSE)
H5 <- "_workspace/gse302339_axisA_20260903/10x_h5"; N_PER <- 300; set.seed(20261018)
cells <- fread("细胞定位/SYN01_Entrez_downstream_20261001/results/lung/objects/axisA_cell_scores.tsv.gz", select = c("gsm", "barcode", "subject_id", "broad_celltype", "passes_QC"))
cells <- cells[passes_QC == TRUE & broad_celltype != "Other/Progenitor"]
ledger <- fread("01_current_analysis/01_storyline_current/06_single_cell_projection/results/GSE302339_axisA_official_ucell_20260903_candidate/tables/sample_ledger.tsv")[matrix_available == TRUE]
# donor-balanced draw: random order within donor, then take cells round-robin across donors
cells[, u := runif(.N)]; setorder(cells, broad_celltype, subject_id, u); cells[, r := seq_len(.N), by = .(broad_celltype, subject_id)]
cells[, v := runif(.N)]; setorder(cells, broad_celltype, r, v); sel <- cells[, head(.SD, N_PER), by = broad_celltype]
fwrite(sel[, .N, by = .(broad_celltype)], file.path(OUT, "qa/lung_reference_cells_per_type.tsv"), sep = "\t")
fwrite(sel[, .(donors = uniqueN(subject_id), max_cells_one_donor = max(table(subject_id))), by = broad_celltype], file.path(OUT, "qa/lung_reference_donor_balance.tsv"), sep = "\t")
read_raw <- function(path) { f <- H5File$new(path, "r"); on.exit(f$close_all())
  m <- new("dgCMatrix", i = as.integer(f[["matrix/indices"]][]), p = as.integer(f[["matrix/indptr"]][]), x = as.numeric(f[["matrix/data"]][]), Dim = as.integer(f[["matrix/shape"]][]))
  rownames(m) <- as.character(f[["matrix/features/id"]][]); colnames(m) <- as.character(f[["matrix/barcodes"]][]); m }
parts <- list(); genes <- NULL
for (k in seq_len(nrow(ledger))) { s <- sel[gsm == ledger$gsm[k]]; if (!nrow(s)) next
  m <- read_raw(file.path(H5, ledger$matrix_file[k])); if (is.null(genes)) genes <- rownames(m); stopifnot(identical(rownames(m), genes))
  j <- match(s$barcode, colnames(m)); stopifnot(!anyNA(j)); x <- m[, j, drop = FALSE]; colnames(x) <- paste(s$gsm, s$barcode, sep = "_"); parts[[k]] <- x }
X <- do.call(cbind, parts); lab <- sel$broad_celltype[match(colnames(X), paste(sel$gsm, sel$barcode, sep = "_"))]
# Ensembl -> Entrez (unique both ways) -> current symbol
mp <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", columns = "ENTREZID"))[!is.na(ENTREZID)]
mp <- mp[!ENSEMBL %in% ENSEMBL[duplicated(ENSEMBL)] & !ENTREZID %in% ENTREZID[duplicated(ENTREZID)]]
mp[, SYMBOL := mapIds(org.Hs.eg.db, ENTREZID, "SYMBOL", "ENTREZID")]; mp <- mp[!is.na(SYMBOL) & !SYMBOL %in% SYMBOL[duplicated(SYMBOL)]]
X <- X[mp$ENSEMBL, ]; rownames(X) <- mp$SYMBOL; keep <- rowSums(X) > 0; X <- X[keep, ]; ent <- mp$ENTREZID[keep]
axisA <- fread("03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisA_confirmed_genes.csv", colClasses = "character")$entrez_id
labs <- gsub("[/ &]+", "_", lab)
write_ref <- function(M, file) { con <- file(file, "w"); writeLines(paste(c("GeneSymbol", labs), collapse = "\t"), con)
  M <- as.matrix(M); storage.mode(M) <- "integer"; write.table(M, con, sep = "\t", quote = FALSE, col.names = FALSE); close(con) }
write_ref(X, file.path(OUT, "upload/lung_scref_A.txt")); write_ref(X[!ent %in% axisA, ], file.path(OUT, "upload/lung_scref_B.txt"))
saveRDS(list(genes_symbol = rownames(X), entrez = ent, labels = labs, cells = colnames(X)), file.path(OUT, "source_data/lung_reference_meta.rds"))
cat("cells", ncol(X), " genes A", nrow(X), " genes B", sum(!ent %in% axisA), " Axis A removed", sum(ent %in% axisA), " datapoints A", ncol(X) * nrow(X), "\n"); print(table(labs))

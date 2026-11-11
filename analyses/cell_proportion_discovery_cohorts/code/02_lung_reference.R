# GSE302339 reference: per broad cell type, mean of per-cell CPM over all QC cells; log2(mean + 1); Ensembl -> Entrez.
.libPaths(c('_workspace/r_libs_deconv_20261002','_workspace/r_libs',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(Matrix);library(hdf5r);library(org.Hs.eg.db)})
OUT <- '03_candidate_analysis/Deconvolution_composition_adjustment_20261002'
H5 <- '_workspace/gse302339_axisA_20260903/10x_h5'
cells <- fread('cell_localization/SYN01_Entrez_downstream_20261001/results/lung/objects/axisA_cell_scores.tsv.gz',select=c('gsm','barcode','broad_celltype','passes_QC'))
stopifnot(all(cells$passes_QC), nrow(cells)==211464)
ledger <- fread('01_current_analysis/01_storyline_current/06_single_cell_projection/results/GSE302339_axisA_official_ucell_20260903_candidate/tables/sample_ledger.tsv')[matrix_available==TRUE]
types <- sort(unique(cells$broad_celltype))
read_raw <- function(path) {
 f <- H5File$new(path,'r'); on.exit(f$close_all())
 m <- new('dgCMatrix',i=as.integer(f[['matrix/indices']][]),p=as.integer(f[['matrix/indptr']][]),x=as.numeric(f[['matrix/data']][]),Dim=as.integer(f[['matrix/shape']][]))
 rownames(m) <- as.character(f[['matrix/features/id']][]); colnames(m) <- as.character(f[['matrix/barcodes']][]); m
}
genes <- NULL; S <- NULL; N <- setNames(numeric(length(types)),types)
for (k in seq_len(nrow(ledger))) {
 meta <- cells[gsm==ledger$gsm[k]]; m <- read_raw(file.path(H5,ledger$matrix_file[k]))
 if (is.null(genes)) { genes <- rownames(m); S <- matrix(0,length(genes),length(types),dimnames=list(genes,types)) }
 stopifnot(identical(rownames(m),genes))
 # per-cell CPM uses the cell's total UMI over all features (as in the SYN01 per-cell QC 'sum')
 tot <- colSums(m); stopifnot(all(tot>0))
 x <- m[,match(meta$barcode,colnames(m))]; stopifnot(!anyNA(match(meta$barcode,colnames(m))))
 cpm <- x %*% Diagonal(x=1e6/tot[meta$barcode])
 G <- sparseMatrix(i=seq_len(nrow(meta)),j=match(meta$broad_celltype,types),x=1,dims=c(nrow(meta),length(types)))
 S <- S + as.matrix(cpm %*% G); N <- N + colSums(G)
 message(k,'/',nrow(ledger))
}
stopifnot(sum(N)==211464)
ref <- log2(sweep(S,2,N,'/')+1)
e2e <- AnnotationDbi::select(org.Hs.eg.db,keys=genes,keytype='ENSEMBL',columns='ENTREZID')
e2e <- e2e[!is.na(e2e$ENTREZID),]
e2e <- e2e[!e2e$ENSEMBL %in% e2e$ENSEMBL[duplicated(e2e$ENSEMBL)] & !e2e$ENTREZID %in% e2e$ENTREZID[duplicated(e2e$ENTREZID)],]
refE <- ref[e2e$ENSEMBL,]; rownames(refE) <- e2e$ENTREZID
saveRDS(list(ref_log2=refE,n_cells=N,genes_ensembl=length(genes),genes_entrez=nrow(refE)),file.path(OUT,'source_data/lung_reference_GSE302339.rds'))
print(N); cat('genes',length(genes),'-> Entrez',nrow(refE),'\n')

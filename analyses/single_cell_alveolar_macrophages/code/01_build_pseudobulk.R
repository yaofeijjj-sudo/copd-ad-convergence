# Sum raw UMI counts per library for three lung cell classes (GSE302339), Ensembl-level.
# Classes: AM (primary), AM family (sensitivity), Myeloid (reference; must equal the 10-14 pseudobulk).
# Reader adapted from AxisA_pseudobulk_disease_20261014/code/01_build_lung_pseudobulk.R (cell filter changed only).
suppressPackageStartupMessages({library(data.table);library(Matrix);library(hdf5r)})
ROOT <- '.'
OUT <- file.path(ROOT,'03_candidate_analysis/AxisA_AM_pseudobulk_20261028')
H5 <- file.path(ROOT,'_workspace/gse302339_axisA_20260903/10x_h5')
CELLS <- file.path(ROOT,'细胞定位/SYN01_Entrez_downstream_20261001/results/lung/objects/axisA_cell_scores.tsv.gz')
LEDGER <- file.path(ROOT,'01_current_analysis/01_storyline_current/06_single_cell_projection/results/GSE302339_axisA_official_ucell_20260903_candidate/tables/sample_ledger.tsv')
REF_MYE <- file.path(ROOT,'03_candidate_analysis/AxisA_pseudobulk_disease_20261014/source_data/lung_myeloid_pseudobulk_counts.tsv.gz')

read_raw <- function(path) {
 f <- H5File$new(path,'r'); on.exit(f$close_all())
 m <- new('dgCMatrix',i=as.integer(f[['matrix/indices']][]),p=as.integer(f[['matrix/indptr']][]),
          x=as.numeric(f[['matrix/data']][]),Dim=as.integer(f[['matrix/shape']][]))
 rownames(m) <- as.character(f[['matrix/features/id']][]); colnames(m) <- as.character(f[['matrix/barcodes']][])
 stopifnot(!anyDuplicated(rownames(m)), !anyDuplicated(colnames(m)))
 m
}

AM_FAMILY <- c('Alveolar macrophages','Alveolar Mph CCL3+','Alveolar Mph proliferating','Alveolar Mph MT-positive')
cells <- fread(CELLS,select=c('gsm','barcode','subject_id','broad_celltype','celltypist_label','passes_QC'))
stopifnot(all(cells$passes_QC))
cells[,`:=`(AM=celltypist_label=='Alveolar macrophages', AMfam=celltypist_label %in% AM_FAMILY, Myeloid=broad_celltype=='Myeloid')]
stopifnot(sum(cells$AM)==64542, all(cells[AMfam==TRUE]$Myeloid))
ledger <- fread(LEDGER)[matrix_available==TRUE]
stopifnot(setequal(unique(cells$gsm),ledger$gsm), nrow(ledger)==65)

classes <- c('AM','AMfam','Myeloid')
genes <- NULL; pb <- list(); ncell <- list()
for (k in seq_len(nrow(ledger))) {
 g <- ledger$gsm[k]; meta <- cells[gsm==g]
 m <- read_raw(file.path(H5,ledger$matrix_file[k]))
 if (is.null(genes)) { genes <- rownames(m); for (cl in classes) pb[[cl]] <- matrix(0,length(genes),nrow(ledger),dimnames=list(genes,ledger$gsm)) }
 stopifnot(identical(rownames(m),genes), uniqueN(meta$subject_id)==1)
 for (cl in classes) {
  bc <- meta[get(cl)==TRUE]$barcode
  if (!length(bc)) next
  idx <- match(bc,colnames(m)); stopifnot(!anyNA(idx))
  x <- m[,idx,drop=FALSE]; stopifnot(all(x@x>=0), all(x@x==floor(x@x)))
  pb[[cl]][,g] <- as.numeric(rowSums(x))
 }
 message(k,'/',nrow(ledger),' ',g,' AM=',sum(meta$AM),' myeloid=',sum(meta$Myeloid))
}
lib <- cells[,.(n_AM=sum(AM),n_AMfam=sum(AMfam),n_Myeloid=sum(Myeloid)),by=.(gsm,subject_id)]
lib <- lib[match(ledger$gsm,gsm)]
fwrite(lib,file.path(OUT,'source_data/library_cell_counts.tsv'),sep='\t')
for (cl in classes) {
 stopifnot(sum(lib[[paste0('n_',cl)]])==sum(cells[[cl]]))
 fwrite(data.table(gene_id=genes,pb[[cl]]),file.path(OUT,sprintf('source_data/library_pseudobulk_%s.tsv.gz',cl)),sep='\t')
}

# Reference check: donor-summed Myeloid must equal the frozen 10-14 pseudobulk exactly.
ref <- fread(REF_MYE); R <- as.matrix(ref[,-1]); rownames(R) <- ref$gene_id
D <- sapply(colnames(R),function(d) rowSums(pb$Myeloid[,lib[subject_id==d]$gsm,drop=FALSE]))
stopifnot(identical(rownames(R),genes), max(abs(D-R))==0)
message('Myeloid pseudobulk identical to 10-14 reference: ',ncol(R),' donors')

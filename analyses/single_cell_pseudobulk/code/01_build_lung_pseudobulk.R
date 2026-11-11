# Sum raw UMI counts of QC-passed Myeloid cells per donor (GSE302339), Ensembl-level.
suppressPackageStartupMessages({library(data.table);library(Matrix);library(hdf5r)})
ROOT <- '.'
OUT <- file.path(ROOT,'03_candidate_analysis/AxisA_pseudobulk_disease_20261014')
H5 <- file.path(ROOT,'_workspace/gse302339_axisA_20260903/10x_h5')
CELLS <- file.path(ROOT,'cell_localization/SYN01_Entrez_downstream_20261001/results/lung/objects/axisA_cell_scores.tsv.gz')
LEDGER <- file.path(ROOT,'01_current_analysis/01_storyline_current/06_single_cell_projection/results/GSE302339_axisA_official_ucell_20260903_candidate/tables/sample_ledger.tsv')
dir.create(file.path(OUT,'source_data'),recursive=TRUE,showWarnings=FALSE)

read_raw <- function(path) {
 f <- H5File$new(path,'r'); on.exit(f$close_all())
 m <- new('dgCMatrix',i=as.integer(f[['matrix/indices']][]),p=as.integer(f[['matrix/indptr']][]),
          x=as.numeric(f[['matrix/data']][]),Dim=as.integer(f[['matrix/shape']][]))
 rownames(m) <- as.character(f[['matrix/features/id']][]); colnames(m) <- as.character(f[['matrix/barcodes']][])
 stopifnot(!anyDuplicated(rownames(m)), !anyDuplicated(colnames(m)))
 m
}

cells <- fread(CELLS,select=c('gsm','barcode','subject_id','disease_group','broad_celltype','passes_QC','age','sex','smoke','packyear'))
stopifnot(all(cells$passes_QC))
mye <- cells[broad_celltype=='Myeloid']
ledger <- fread(LEDGER)[matrix_available==TRUE]
stopifnot(setequal(unique(cells$gsm),ledger$gsm))

donors <- sort(unique(mye$subject_id))
genes <- NULL; pb <- NULL; n_cells <- setNames(integer(length(donors)),donors); n_lib <- n_cells
for (k in seq_len(nrow(ledger))) {
 g <- ledger$gsm[k]; meta <- mye[gsm==g]
 if (!nrow(meta)) next
 m <- read_raw(file.path(H5,ledger$matrix_file[k]))
 if (is.null(genes)) { genes <- rownames(m); pb <- matrix(0,length(genes),length(donors),dimnames=list(genes,donors)) }
 stopifnot(identical(rownames(m),genes))
 idx <- match(meta$barcode,colnames(m)); stopifnot(!anyNA(idx))
 x <- m[,idx,drop=FALSE]
 stopifnot(all(x@x>=0), all(x@x==floor(x@x)))
 d <- unique(meta$subject_id); stopifnot(length(d)==1)  # one donor per library
 pb[,d] <- pb[,d] + as.numeric(rowSums(x))
 n_cells[d] <- n_cells[d] + ncol(x); n_lib[d] <- n_lib[d] + 1L
 message(k,'/',nrow(ledger),' ',g,' myeloid=',ncol(x))
}
stopifnot(sum(n_cells)==nrow(mye), all(n_cells>=20), length(donors)==56)
donor_meta <- unique(mye[,.(subject_id,disease_group,age,sex,smoke,packyear)])
stopifnot(!anyDuplicated(donor_meta$subject_id))
donor_meta[,`:=`(n_cells=n_cells[subject_id],n_libraries=n_lib[subject_id])]
setorder(donor_meta,subject_id)
fwrite(donor_meta,file.path(OUT,'source_data/lung_myeloid_donors.tsv'),sep='\t')
dt <- data.table(gene_id=genes,pb); fwrite(dt,file.path(OUT,'source_data/lung_myeloid_pseudobulk_counts.tsv.gz'),sep='\t')
message('done: ',sum(n_cells),' cells, ',length(donors),' donors, ',length(genes),' genes')

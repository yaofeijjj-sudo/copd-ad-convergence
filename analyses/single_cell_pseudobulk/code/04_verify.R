# QA: (1) lung donor library sizes vs independent per-cell QC totals; (2) fry vs exact-rotation roast;
# (3) input identities and counts recorded in qa/verification.json.
suppressPackageStartupMessages({library(data.table);library(edgeR);library(jsonlite)})
ROOT <- '.'
OUT <- file.path(ROOT,'03_candidate_analysis/AxisA_pseudobulk_disease_20261014')
SYN <- file.path(ROOT,'细胞定位/SYN01_Entrez_downstream_20261001/results')
dir.create(file.path(OUT,'qa'),showWarnings=FALSE)
sha <- function(f) strsplit(system2('shasum',c('-a','256',shQuote(f)),stdout=TRUE),' ')[[1]][1]
z <- function(x) (x-mean(x))/sd(x)

cells <- fread(file.path(SYN,'lung/objects/axisA_cell_scores.tsv.gz'),select=c('subject_id','broad_celltype','sum'))
qc_tot <- cells[broad_celltype=='Myeloid',.(qc_sum=sum(sum),n=.N),by=subject_id]
lc <- fread(file.path(OUT,'source_data/lung_myeloid_pseudobulk_counts.tsv.gz'))
pb_tot <- data.table(subject_id=names(lc)[-1],pb_sum=colSums(as.matrix(lc[,-1])))
chk <- merge(qc_tot,pb_tot,by='subject_id')
lung_total_maxdiff <- max(abs(chk$qc_sum-chk$pb_sum))

# fry vs roast (lung primary, brain primary), same filtered/normalised objects as 03
set.seed(20261014)
lm <- fread(file.path(OUT,'source_data/lung_myeloid_donors.tsv'))
L <- as.matrix(lc[,-1]); rownames(L) <- lc$gene_id; L <- L[,lm$subject_id]
lm[,`:=`(disease=factor(disease_group,levels=c('Control','COPD')),age_z=z(age),sex=factor(sex),smoke=factor(smoke,levels=c('Never','Former','Current')))]
dl <- model.matrix(~disease+age_z+sex+smoke,data=lm)
bc <- fread(file.path(ROOT,'虚拟敲除/AxisA_疾病标尺核对_20260925/source_data/brain_donor_pseudobulk_counts.tsv.gz'))
map <- fread(file.path(OUT,'source_data/brain_symbol_to_ensembl.tsv'))
bm <- fread(file.path(SYN,'disease/brain_CPS_donor_data.tsv'))[Supertype=='Micro-PVM_2'];setorder(bm,`Donor ID`)
B <- as.matrix(bc[,-1]); rownames(B) <- map$gene_id; B <- B[,bm$`Donor ID`]
bm[,`:=`(CPS_Local_z=z(cps_local),age_z=z(age),PMI_z=z(pmi),RIN_z=z(rin),sex=factor(sex),study=factor(study))]
db <- model.matrix(~CPS_Local_z+age_z+sex+PMI_z+RIN_z+study,data=bm)
cmp <- function(counts,design,coef,axis) {
 y <- DGEList(counts); y <- normLibSizes(y[filterByExpr(y,design=design),,keep.lib.sizes=FALSE]); y <- estimateDisp(y,design)
 idx <- list(AxisA=intersect(axis,rownames(y))); j <- match(coef,colnames(design))
 f <- fry(y,index=idx,design=design,contrast=j); r <- roast(y,index=idx,design=design,contrast=j,nrot=99999)
 data.table(fry_direction=f$Direction,fry_p=f$PValue,roast_direction=r$Direction,roast_p=r$PValue)
}
roast_chk <- rbind(cbind(organ='lung',cmp(L,dl,'diseaseCOPD',readLines(file.path(SYN,'lung_axis_measured.txt')))),
                   cbind(organ='brain',cmp(B,db,'CPS_Local_z',readLines(file.path(SYN,'brain_axis_measured.txt')))))
fwrite(roast_chk,file.path(OUT,'qa/fry_vs_roast.tsv'),sep='\t')
print(roast_chk)

qa <- list(
 lung=list(cells=sum(lm$n_cells),donors=nrow(lm),COPD=sum(lm$disease_group=='COPD'),Control=sum(lm$disease_group=='Control'),
           matches_S3_Table=(sum(lm$n_cells)==116558 && nrow(lm)==56),
           donor_total_UMI_vs_cell_QC_sum_max_abs_diff=lung_total_maxdiff,
           counts_sha256=sha(file.path(OUT,'source_data/lung_myeloid_pseudobulk_counts.tsv.gz'))),
 brain=list(nuclei=sum(bm$n_nuclei),donors=nrow(bm),matches_S3_Table=(sum(bm$n_nuclei)==23198 && nrow(bm)==79),
           counts_sha256=sha(file.path(ROOT,'虚拟敲除/AxisA_疾病标尺核对_20260925/source_data/brain_donor_pseudobulk_counts.tsv.gz')),
           expected_sha256_from_20260925='b997c1f984e9080cb676f1b0c09bb365e9b124ea22a318c78aaa040a15236b5b'),
 axis_lists=list(lung=sha(file.path(SYN,'lung_axis_measured.txt')),brain=sha(file.path(SYN,'brain_axis_measured.txt'))),
 fry_vs_roast=roast_chk,
 software=list(R=R.version.string,edgeR=as.character(packageVersion('edgeR')),limma=as.character(packageVersion('limma'))))
writeLines(toJSON(qa,auto_unbox=TRUE,pretty=TRUE,digits=NA),file.path(OUT,'qa/verification.json'))
cat('lung total diff', lung_total_maxdiff,'\n')

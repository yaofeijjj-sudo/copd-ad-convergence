# Donor-level COPD vs control tests of Axis A within alveolar macrophages.
# Pseudobulk: edgeR filterByExpr + TMM + estimateDisp + fry; voom + CAMERA (inter.gene.cor = NA).
# UCell: library median -> donor median; two-sided Wilcoxon with Hodges-Lehmann shift and 95% CI.
suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
ROOT <- '.'
OUT <- file.path(ROOT,'03_candidate_analysis/AxisA_AM_pseudobulk_20261028')
SYN <- file.path(ROOT,'cell_localization/SYN01_Entrez_downstream_20261001/results')
CELLS <- file.path(SYN,'lung/objects/axisA_cell_scores.tsv.gz')
FLAG <- file.path(ROOT,'03_candidate_analysis/Lung_sc_library_QC_check_20261028/results/library_qc_with_flag.tsv')
z <- function(x) (x-mean(x))/sd(x)
AM_FAMILY <- c('Alveolar macrophages','Alveolar Mph CCL3+','Alveolar Mph proliferating','Alveolar Mph MT-positive')

cells <- fread(CELLS,select=c('gsm','subject_id','disease_group','broad_celltype','celltypist_label','AxisA_UCell','age','sex','smoke'))
cells[,`:=`(AM=celltypist_label=='Alveolar macrophages', AMfam=celltypist_label %in% AM_FAMILY, Myeloid=broad_celltype=='Myeloid')]
donors <- unique(cells[,.(subject_id,disease_group,age,sex,smoke)]); stopifnot(!anyDuplicated(donors$subject_id))
donors[,`:=`(disease=factor(disease_group,levels=c('Control','COPD')),sex=factor(sex),smoke=factor(smoke,levels=c('Never','Former','Current')))]
lib <- fread(file.path(OUT,'source_data/library_cell_counts.tsv'))
flagged <- fread(FLAG)[flagged==TRUE]$gsm; stopifnot(setequal(flagged,c('GSM9102301','GSM9102263')))
axis <- readLines(file.path(SYN,'lung_axis_measured.txt')); stopifnot(length(axis)==136)

donor_pb <- function(cl, libs) {
 x <- fread(file.path(OUT,sprintf('source_data/library_pseudobulk_%s.tsv.gz',cl)))
 M <- as.matrix(x[,-1]); rownames(M) <- x$gene_id
 l <- lib[gsm %in% libs]
 ds <- sort(unique(l$subject_id))
 P <- sapply(ds,function(d) rowSums(M[,l[subject_id==d]$gsm,drop=FALSE]))
 n <- sapply(ds,function(d) sum(l[subject_id==d][[paste0('n_',cl)]]))
 list(counts=P,n=n)
}

run_pb <- function(label, cl, libs=lib$gsm, min_cells=20) {
 d <- donor_pb(cl, libs); keep_d <- names(d$n)[d$n>=min_cells]
 C <- d$counts[,keep_d,drop=FALSE]; stopifnot(all(C==round(C)), all(C>=0))
 meta <- donors[match(keep_d,subject_id)]; meta[,age_z:=z(age)]
 stopifnot(!anyNA(meta[,.(disease,age_z,sex,smoke)]))
 des <- model.matrix(~disease+age_z+sex+smoke,data=meta); rownames(des) <- keep_d
 stopifnot(qr(des)$rank==ncol(des))
 y <- DGEList(C); keep <- filterByExpr(y,design=des); y <- normLibSizes(y[keep,,keep.lib.sizes=FALSE])
 members <- intersect(axis,rownames(y)); j <- match('diseaseCOPD',colnames(des))
 y <- estimateDisp(y,des); fr <- fry(y,index=list(AxisA=members),design=des,contrast=j)
 v <- voom(y,des); cm <- camera(v,index=list(AxisA=which(rownames(v) %in% members)),design=des,contrast=j,inter.gene.cor=NA)
 tt <- topTags(glmQLFTest(glmQLFit(y,des),coef=j),n=Inf,sort.by='none')$table
 g <- data.table(analysis=label,gene_id=rownames(tt),tt[,c('logFC','logCPM','PValue','FDR')])[gene_id %in% members]
 row <- data.table(analysis=label,class=cl,n_donors=ncol(C),n_COPD=sum(meta$disease=='COPD'),n_control=sum(meta$disease=='Control'),
   n_cells=sum(d$n[keep_d]),axis_after_filter=length(members),fry_direction=fr$Direction,fry_p=fr$PValue,
   camera_direction=cm$Direction,camera_p=cm$PValue,camera_corr=cm$Correlation,
   members_up=sum(g$logFC>0),members_down=sum(g$logFC<0),median_member_logFC=median(g$logFC))
 list(row=row,genes=g)
}

run_ucell <- function(label, cl, libs=lib$gsm, min_cells=20) {
 lm <- cells[get(cl)==TRUE & gsm %in% libs,.(n=.N,med=median(AxisA_UCell)),by=.(gsm,subject_id,disease_group)]
 dm <- lm[,.(n=sum(n),med=median(med)),by=.(subject_id,disease_group)][n>=min_cells]
 x <- dm[disease_group=='COPD']$med; y <- dm[disease_group=='Control']$med
 w <- suppressWarnings(wilcox.test(x,y,exact=FALSE,conf.int=TRUE))
 data.table(analysis=label,class=cl,n_donors=nrow(dm),n_COPD=length(x),n_control=length(y),
   median_control=median(y),median_COPD=median(x),diff_medians=median(x)-median(y),
   HL_shift=unname(w$estimate),HL_low=w$conf.int[1],HL_high=w$conf.int[2],wilcoxon_p=w$p.value)
}

unflag <- setdiff(lib$gsm,flagged)
PB <- list(
 run_pb('primary: AM','AM'),
 run_pb('sensitivity a: AM, >=50 cells','AM',min_cells=50),
 run_pb('sensitivity b: AM, flagged libraries removed','AM',libs=unflag),
 run_pb('sensitivity c: AM family','AMfam'),
 run_pb('reference: Myeloid (S6 Table)','Myeloid'))
UC <- rbindlist(list(
 run_ucell('primary: AM','AM'),
 run_ucell('sensitivity a: AM, >=50 cells','AM',min_cells=50),
 run_ucell('sensitivity b: AM, flagged libraries removed','AM',libs=unflag),
 run_ucell('sensitivity c: AM family','AMfam'),
 run_ucell('reference: Myeloid (S6 Table)','Myeloid')))
pbr <- rbindlist(lapply(PB,`[[`,'row')); genes <- rbindlist(lapply(PB,`[[`,'genes'))

# Reference rows must reproduce S6 Table.
ref_pb <- pbr[class=='Myeloid']; ref_uc <- UC[class=='Myeloid']
chk <- c(fry=abs(ref_pb$fry_p-0.1151250275), camera=abs(ref_pb$camera_p-0.0561733203), wilcoxon=abs(ref_uc$wilcoxon_p-0.01032))
print(chk); stopifnot(chk[1:2]<1e-6, chk[3]<5e-5)

# Descriptive: AM share of myeloid cells per donor (not tested).
share <- lib[,.(n_AM=sum(n_AM),n_Myeloid=sum(n_Myeloid)),by=subject_id][donors,on='subject_id']
share[,AM_share:=n_AM/n_Myeloid]
share_sum <- share[,.(n_donors=.N,median_AM_share=median(AM_share),q25=quantile(AM_share,.25),q75=quantile(AM_share,.75)),by=disease_group]

fwrite(pbr,file.path(OUT,'results/pseudobulk_fry_camera.tsv'),sep='\t')
fwrite(UC,file.path(OUT,'results/donor_ucell_wilcoxon.tsv'),sep='\t')
fwrite(genes,file.path(OUT,'results/axisA_member_logFC.tsv'),sep='\t')
fwrite(share[,.(subject_id,disease_group,n_AM,n_Myeloid,AM_share)],file.path(OUT,'results/AM_share_of_myeloid_by_donor.tsv'),sep='\t')
fwrite(share_sum,file.path(OUT,'results/AM_share_of_myeloid_summary.tsv'),sep='\t')
print(pbr[,.(analysis,n_donors,axis_after_filter,fry_direction,fry_p=signif(fry_p,3),camera_direction,camera_p=signif(camera_p,3),camera_corr=round(camera_corr,3),members_up,members_down)])
print(UC[,.(analysis,n_donors,median_control=signif(median_control,4),median_COPD=signif(median_COPD,4),HL_shift=signif(HL_shift,3),HL_low=signif(HL_low,3),HL_high=signif(HL_high,3),wilcoxon_p=signif(wilcoxon_p,3))])
print(share_sum)
writeLines(capture.output(sessionInfo()),file.path(OUT,'logs/R_sessionInfo.txt'))

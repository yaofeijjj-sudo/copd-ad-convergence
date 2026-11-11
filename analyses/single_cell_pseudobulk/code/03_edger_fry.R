# Donor pseudobulk: edgeR filterByExpr + TMM + fry (self-contained, directional) for the fixed Axis A set.
# Follows DESIGN.md; edgeR v4 official workflow (Chen et al. 2025 NAR).
suppressPackageStartupMessages({library(data.table);library(edgeR)})
ROOT <- '.'
OUT <- file.path(ROOT,'03_candidate_analysis/AxisA_pseudobulk_disease_20261014')
SYN <- file.path(ROOT,'细胞定位/SYN01_Entrez_downstream_20261001/results')
dir.create(file.path(OUT,'results'),showWarnings=FALSE)
z <- function(x) (x-mean(x))/sd(x)

run_one <- function(organ,model,counts,samples,design,coef,axis) {
 stopifnot(identical(colnames(counts),rownames(design)), qr(design)$rank==ncol(design))
 y <- DGEList(counts)
 keep <- filterByExpr(y,design=design)
 y <- normLibSizes(y[keep,,keep.lib.sizes=FALSE])
 members <- intersect(axis,rownames(y))
 j <- match(coef,colnames(design))
 y <- estimateDisp(y,design)  # fry.DGEList requires dispersions (edgeR z-score transformation)
 fr <- fry(y,index=list(AxisA=members),design=design,contrast=j)
 fit <- glmQLFit(y,design)
 tt <- topTags(glmQLFTest(fit,coef=j),n=Inf,sort.by='none')$table
 tt$gene_id <- rownames(tt); tt$axisA <- tt$gene_id %in% members
 genes <- data.table(organ=organ,model=model,tt[tt$axisA,c('gene_id','logFC','logCPM','PValue','FDR')])
 row <- data.table(organ=organ,model=model,coefficient=coef,n_donors=ncol(counts),
   n_genes_tested=nrow(y),axis_measured=length(axis),axis_after_filter=length(members),
   fry_direction=fr$Direction,fry_p=fr$PValue,fry_p_mixed=fr$PValue.Mixed,
   members_up=sum(genes$logFC>0),members_down=sum(genes$logFC<0),
   median_member_logFC=median(genes$logFC),
   members_genomewide_FDR_lt_0.05=sum(genes$FDR<0.05))
 list(row=row,genes=genes)
}

## Lung: GSE302339 Myeloid, COPD vs control
lc <- fread(file.path(OUT,'source_data/lung_myeloid_pseudobulk_counts.tsv.gz'))
lm <- fread(file.path(OUT,'source_data/lung_myeloid_donors.tsv'))
L <- as.matrix(lc[,-1]); rownames(L) <- lc$gene_id; L <- L[,lm$subject_id]
stopifnot(all(L==round(L)), all(L>=0))
lm[,`:=`(disease=factor(disease_group,levels=c('Control','COPD')),age_z=z(age),
         sex=factor(sex),smoke=factor(smoke,levels=c('Never','Former','Current')))]
stopifnot(!anyNA(lm[,.(disease,age_z,sex,smoke)]))
lung_axis <- readLines(file.path(SYN,'lung_axis_measured.txt'))
dl1 <- model.matrix(~disease+age_z+sex+smoke,data=lm); rownames(dl1) <- lm$subject_id
dl2 <- model.matrix(~disease,data=lm); rownames(dl2) <- lm$subject_id
r1 <- run_one('lung','primary: ~ disease + age + sex + smoking status',L,lm,dl1,'diseaseCOPD',lung_axis)
r2 <- run_one('lung','secondary: ~ disease',L,lm,dl2,'diseaseCOPD',lung_axis)

## Brain: SEA-AD DFC Micro-PVM_2, per 1 SD CPS_Local (Figure 5 model)
bc <- fread(file.path(ROOT,'虚拟敲除/AxisA_疾病标尺核对_20260925/source_data/brain_donor_pseudobulk_counts.tsv.gz'))
map <- fread(file.path(OUT,'source_data/brain_symbol_to_ensembl.tsv'))
stopifnot(identical(bc$gene,map$gene))
bm <- fread(file.path(SYN,'disease/brain_CPS_donor_data.tsv'))[Supertype=='Micro-PVM_2']
stopifnot(nrow(bm)==79, !anyDuplicated(bm$`Donor ID`))
setorder(bm,`Donor ID`)
B <- as.matrix(bc[,-1]); rownames(B) <- map$gene_id; B <- B[,bm$`Donor ID`]
stopifnot(all(B==round(B)), all(B>=0))
bm[,`:=`(CPS_Local_z=z(cps_local),age_z=z(age),PMI_z=z(pmi),RIN_z=z(rin),sex=factor(sex),study=factor(study))]
brain_axis <- readLines(file.path(SYN,'brain_axis_measured.txt'))
db <- model.matrix(~CPS_Local_z+age_z+sex+PMI_z+RIN_z+study,data=bm); rownames(db) <- bm$`Donor ID`
r3 <- run_one('brain','primary: ~ CPS_Local_z + age + sex + PMI + RIN + study',B,bm,db,'CPS_Local_z',brain_axis)

res <- rbindlist(list(r1$row,r2$row,r3$row))
prim <- res$model %like% '^primary'
res[prim,BH_across_two_organs_reference:=p.adjust(fry_p,'BH')]
fwrite(res,file.path(OUT,'results/fry_summary.tsv'),sep='\t')
mg <- rbindlist(list(r1$genes,r2$genes,r3$genes)); mg[,symbol:=map$gene[match(gene_id,map$gene_id)]]; setcolorder(mg,c('organ','model','gene_id','symbol'))
fwrite(mg,file.path(OUT,'results/axisA_member_logFC.tsv'),sep='\t')
fwrite(lm[,.(subject_id,disease_group,age,sex,smoke,packyear,n_cells,n_libraries)],file.path(OUT,'results/lung_donors_used.tsv'),sep='\t')
fwrite(bm[,.(`Donor ID`,cps_local,age,sex,pmi,rin,study,n_nuclei)],file.path(OUT,'results/brain_donors_used.tsv'),sep='\t')
print(res)
writeLines(capture.output(sessionInfo()),file.path(OUT,'logs/R_sessionInfo.txt'))

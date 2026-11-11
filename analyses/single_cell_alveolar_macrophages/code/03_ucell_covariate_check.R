# Pre-specified rule: UCell (Wilcoxon P < 0.05) and pseudobulk (P >= 0.05) disagree.
# Check whether the covariates used in the pseudobulk model explain the gap: linear model on donor UCell medians.
suppressPackageStartupMessages(library(data.table))
ROOT <- '.'; OUT <- file.path(ROOT,'03_candidate_analysis/AxisA_AM_pseudobulk_20261028')
cells <- fread(file.path(ROOT,'cell_localization/SYN01_Entrez_downstream_20261001/results/lung/objects/axisA_cell_scores.tsv.gz'),
  select=c('gsm','subject_id','disease_group','broad_celltype','celltypist_label','AxisA_UCell','age','sex','smoke'))
z <- function(x) (x-mean(x))/sd(x)
one <- function(label, sel) {
 lm1 <- cells[sel,.(med=median(AxisA_UCell),n=.N),by=.(gsm,subject_id,disease_group,age,sex,smoke)]
 d <- lm1[,.(med=median(med),n=sum(n)),by=.(subject_id,disease_group,age,sex,smoke)][n>=20]
 d[,`:=`(disease=factor(disease_group,levels=c('Control','COPD')),age_z=z(age),smoke=factor(smoke,levels=c('Never','Former','Current')),y=z(med))]
 rbindlist(lapply(list(unadjusted=y~disease, adjusted=y~disease+age_z+sex+smoke), function(f) {
  co <- summary(lm(f,data=d))$coefficients['diseaseCOPD',]; ci <- confint(lm(f,data=d))['diseaseCOPD',]
  data.table(class=label,model=deparse(f),n_donors=nrow(d),effect_SD=co[1],ci_low=ci[1],ci_high=ci[2],p=co[4])}))
}
res <- rbind(one('AM',cells$celltypist_label=='Alveolar macrophages'), one('Myeloid',cells$broad_celltype=='Myeloid'))
fwrite(res,file.path(OUT,'results/donor_ucell_lm_covariate_check.tsv'),sep='\t'); print(res)

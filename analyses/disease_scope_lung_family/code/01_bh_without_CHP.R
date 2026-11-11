# Lung comparator family without chronic hypersensitivity pneumonitis (CHP).
# CHP and control were sequenced in non-overlapping processing batches in GSE150910 (disease fully confounded with batch),
# so per purchased advice 3.3 the comparison becomes descriptive and leaves the BH family. Effects and P values are unchanged.
suppressPackageStartupMessages(library(data.table))
ROOT <- '.'
ds <- fread(file.path(ROOT,'05_github_release/copd-ad-convergence/validation/reference/disease_scope.tsv'))
lung <- ds[id %in% c('GSE150910 IPF','GSE150910 smoking')]
stopifnot(nrow(lung)==2)
lung[,`:=`(fdr_with_CHP=fdr, fdr=p.adjust(p,'BH'))]
chp <- ds[id=='GSE150910 CHP'][,`:=`(fdr_with_CHP=fdr, fdr=NA_real_)]
out <- rbind(chp,lung)[,.(id,n_case,n_control,effect_SD,ci_low,ci_high,p,fdr_with_CHP,fdr_without_CHP=fdr)]
fwrite(out,file.path(ROOT,'03_candidate_analysis/AxisA_disease_scope_lung_family_20261028/results/lung_family_without_CHP.tsv'),sep='\t',na='NA')
print(out)

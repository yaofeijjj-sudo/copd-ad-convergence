# Descriptive QA: correlations among estimated proportions and VIF of the multi-type models (no new inference).
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
OUT <- '03_candidate_analysis/Deconvolution_composition_adjustment_20261002'
p <- readRDS(file.path(OUT,'results/estimated_proportions.rds'))
e1 <- new.env(); load('WGCNA_figures/reproduction_materials/analysis/01_lung_wgcna/results/COPD_WGCNA/copd_wgcna_final_env.RData',envir=e1)
e2 <- new.env(); load('01_Results_and_Methods/01_Figure1_WGCNA/analysis/02_brain_wgcna/results/AD_WGCNA_no_residual_after_age_outlier/ad_wgcna_final_env.RData',envir=e2)
vif <- function(X) { X <- as.data.frame(scale(X)); sapply(names(X), function(v) { r2 <- summary(lm(X[[v]] ~ ., data=X[setdiff(names(X),v)]))$r.squared; 1/(1-r2) }) }
lb <- p$lung$B; bb <- p$brain$B
lungX <- cbind(disease=e1$clinical_train$copd_status, e1$clinical_train[,c('age_z','sex_male','packyears_z')], lb[,setdiff(colnames(lb),names(which.max(colMeans(lb))))])
brainX <- cbind(disease=e2$clinical_train$Disease_State, e2$clinical_train[,c('age_z','sex_male')], bb)
out <- list(brain_proportion_correlations=round(cor(bb),3), lung_proportion_correlations=round(cor(lb),3),
            vif_brain_multi_type=round(vif(brainX),2), vif_lung_multi_type=round(vif(lungX),2),
            brain_R2_disease_on_six_surrogates=summary(lm(e2$clinical_train$Disease_State ~ bb))$r.squared,
            lung_R2_disease_on_proportions=summary(lm(e1$clinical_train$copd_status ~ lb))$r.squared)
writeLines(toJSON(out,pretty=TRUE,auto_unbox=TRUE,digits=NA),file.path(OUT,'qa/collinearity_versionB.json'))
print(out)

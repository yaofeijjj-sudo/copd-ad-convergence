# Estimate cell proportions (lung dtangle; brain BRETIGEA), versions A/B/C, then add them to the Figure 2d models.
.libPaths(c('_workspace/r_libs_deconv_20261002','_workspace/r_libs',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(dtangle);library(BRETIGEA);library(org.Hs.eg.db)})
OUT <- '03_candidate_analysis/Deconvolution_composition_adjustment_20261002'
ID <- '03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results'
bulk <- readRDS(file.path(OUT,'source_data/bulk_full_matrices.rds'))
lref <- readRDS(file.path(OUT,'source_data/lung_reference_GSE302339.rds'))
axisA <- fread(file.path(ID,'cross/AxisA_confirmed_genes.csv'),colClasses='character')$entrez_id; stopifnot(length(axisA)==137)
lung_mod <- fread(file.path(ID,'identity/lung_confirmed_module_gene_ids.csv'),colClasses='character')
brain_mod <- fread(file.path(ID,'identity/brain_confirmed_module_gene_ids.csv'),colClasses='character')
red <- unique(lung_mod[module=='red',entrez_id]); blue <- unique(brain_mod[module=='blue',entrez_id])
cat('red',length(red),'blue',length(blue),'\n')
drop_sets <- list(A=character(0),B=axisA)

## ---------- Lung: dtangle ----------
six <- c('B/Plasma','Endothelial','Epithelial','Myeloid','Stromal','T/NK')
R <- lref$ref_log2[,six]
common <- intersect(rownames(bulk$lung),rownames(R))
lung_props <- list(); lung_markers <- list()
for (v in c('A','B','C')) {
 drop <- switch(v,A=character(0),B=axisA,C=red)
 g <- setdiff(common,drop)
 dt <- dtangle(Y=t(bulk$lung[g,]),references=t(R[g,]),data_type='microarray-gene')
 est <- dt$estimates; colnames(est) <- six
 lung_props[[v]] <- est
 lung_markers[[v]] <- data.table(version=v,cell=six,n_markers=lengths(dt$markers),
                                 n_markers_in_AxisA=sapply(dt$markers,function(i) sum(g[i] %in% axisA)),
                                 n_genes_available=length(g))
}

## ---------- Brain: BRETIGEA ----------
e2s <- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(bulk$brain),keytype='ENTREZID',columns='SYMBOL')
e2s <- e2s[!is.na(e2s$SYMBOL) & !duplicated(e2s$ENTREZID) & !e2s$SYMBOL %in% e2s$SYMBOL[duplicated(e2s$SYMBOL)],]
Bs <- bulk$brain[e2s$ENTREZID,]; rownames(Bs) <- e2s$SYMBOL
mk <- as.data.table(markers_df_human_brain)
s2e <- AnnotationDbi::select(org.Hs.eg.db,keys=unique(mk$markers),keytype='SYMBOL',columns='ENTREZID')
amb <- unique(s2e$SYMBOL[duplicated(s2e$SYMBOL)]); s2e <- s2e[!s2e$SYMBOL %in% amb,]
mk[,entrez:=s2e$ENTREZID[match(markers,s2e$SYMBOL)]]
brain_props <- list(); brain_markers <- list()
for (v in c('A','B','C')) {
 drop <- switch(v,A=character(0),B=axisA,C=blue)
 m <- as.data.frame(mk[is.na(entrez) | !entrez %in% drop,.(markers,cell)])
 invisible(capture.output(est <- findCells(inputMat=Bs,markers=m,nMarker=50,method='SVD',scale=TRUE)))
 brain_props[[v]] <- est
 used <- rbindlist(lapply(unique(m$cell),function(ct){ mm <- m[m$cell==ct,]; pres <- mm$markers[mm$markers %in% rownames(Bs)]
   pres <- pres[vapply(pres,function(x) sum(Bs[x,]>0)>0,logical(1))]; data.table(cell=ct,marker=head(pres,50))}))
 used[,entrez:=mk$entrez[match(marker,mk$markers)]]
 brain_markers[[v]] <- used[,.(n_markers=.N,n_markers_in_AxisA=sum(entrez %in% axisA,na.rm=TRUE),
                               n_markers_in_blue=sum(entrez %in% blue,na.rm=TRUE)),by=cell][,version:=v][]
}

## ---------- Models ----------
e1 <- new.env(); load('WGCNA_figures/reproduction_materials/analysis/01_lung_wgcna/results/COPD_WGCNA/copd_wgcna_final_env.RData',envir=e1)
e2 <- new.env(); load('01_Results_and_Methods/01_Figure1_WGCNA/analysis/02_brain_wgcna/results/AD_WGCNA_no_residual_after_age_outlier/ad_wgcna_final_env.RData',envir=e2)
ld <- as.data.table(e1$clinical_train); ld[,`:=`(red=e1$MEs[sample_id,'MEred'],turquoise=e1$MEs[sample_id,'MEturquoise'])]
bd <- as.data.table(e2$clinical_train); bd[,`:=`(blue=e2$MEs[sample_id,'MEblue'],turquoise=e2$MEs[sample_id,'MEturquoise'])]
stopifnot(!anyNA(ld$red),!anyNA(bd$blue),identical(ld$sample_id,colnames(bulk$lung)),identical(bd$sample_id,colnames(bulk$brain)))
fit <- function(d,y,disease,covs,adj=NULL) {
 f <- as.formula(paste0('scale(',y,') ~ ',paste(c(disease,covs,adj),collapse=' + ')))
 s <- summary(lm(f,data=d))$coefficients; ci <- confint(lm(f,data=d))[disease,]
 data.table(beta=s[disease,1],se=s[disease,2],lo=ci[1],hi=ci[2],p=s[disease,4])
}
L_cov <- c('age_z','sex_male','packyears_z'); B_cov <- c('age_z','sex_male')
base <- rbind(cbind(organ='lung',module='red',fit(ld,'red','copd_status',L_cov)),
              cbind(organ='brain',module='blue',fit(bd,'blue','Disease_State',B_cov)),
              cbind(organ='lung',module='turquoise',fit(ld,'turquoise','copd_status',L_cov)),
              cbind(organ='brain',module='turquoise',fit(bd,'turquoise','Disease_State',B_cov)))
fig2d <- fread('cross WGCNA/forest_revision_20260917/source_data/panel_d_adjusted_effects.csv')
base_check <- max(abs(base$beta - fig2d$effect[match(paste(base$organ,base$module),paste(tolower(fig2d$organ),fig2d$module))]))
cat('baseline vs Figure 2d max|diff|:',base_check,'\n'); stopifnot(base_check<1e-6)

res <- list(); propdiff <- list(); cors <- list()
add <- function(organ,module,version,method,adjusted_for,d,y,disease,covs,props,cols,role) {
 for (cc in cols) d[[paste0('p_',make.names(cc))]] <- as.numeric(scale(props[,cc]))
 adj <- paste0('p_',make.names(cols))
 o <- organ; mo <- module; b <- base[organ==o & module==mo]; stopifnot(nrow(b)==1)
 a <- fit(d,y,disease,covs,adj)
 res[[length(res)+1]] <<- data.table(organ,module,role,method,version,adjusted_for=paste(cols,collapse='+'),
   beta_base=b$beta,lo_base=b$lo,hi_base=b$hi,p_base=b$p,beta_adj=a$beta,lo_adj=a$lo,hi_adj=a$hi,p_adj=a$p,
   pct_change=1-a$beta/b$beta, ci_excludes_0=(a$lo>0|a$hi<0), same_direction=sign(a$beta)==sign(b$beta))
}
for (v in c('A','B','C')) {
 role <- if (v=='B') 'primary' else 'sensitivity'
 add('lung','red',v,'dtangle','Myeloid',copy(ld),'red','copd_status',L_cov,lung_props[[v]],'Myeloid',role)
 add('brain','blue',v,'BRETIGEA','mic',copy(bd),'blue','Disease_State',B_cov,brain_props[[v]],'mic',role)
 add('lung','turquoise',v,'dtangle','Myeloid',copy(ld),'turquoise','copd_status',L_cov,lung_props[[v]],'Myeloid',if(v=='B')'secondary' else 'sensitivity')
 add('brain','turquoise',v,'BRETIGEA','neu',copy(bd),'turquoise','Disease_State',B_cov,brain_props[[v]],'neu',if(v=='B')'secondary' else 'sensitivity')
 refct <- names(which.max(colMeans(lung_props[[v]])))
 add('lung','red',v,'dtangle','all types except reference',copy(ld),'red','copd_status',L_cov,lung_props[[v]],setdiff(six,refct),'sensitivity (multi-type)')
 add('brain','blue',v,'BRETIGEA','all six',copy(bd),'blue','Disease_State',B_cov,brain_props[[v]],colnames(brain_props[[v]]),'sensitivity (multi-type)')
 for (cc in six) { d <- copy(ld); d$prop <- as.numeric(scale(lung_props[[v]][,cc]))
   s <- summary(lm(as.formula(paste('prop ~ copd_status +',paste(L_cov,collapse='+'))),data=d))$coefficients
   propdiff[[length(propdiff)+1]] <- data.table(organ='lung',version=v,cell=cc,beta_disease_SD=s['copd_status',1],p=s['copd_status',4],
     mean_prop_control=mean(lung_props[[v]][ld$copd_status==0,cc]),mean_prop_case=mean(lung_props[[v]][ld$copd_status==1,cc]))
   cors[[length(cors)+1]] <- data.table(organ='lung',version=v,cell=cc,r_with_red=cor(ld$red,lung_props[[v]][,cc]),r_with_turquoise=cor(ld$turquoise,lung_props[[v]][,cc])) }
 for (cc in colnames(brain_props[[v]])) { d <- copy(bd); d$prop <- as.numeric(scale(brain_props[[v]][,cc]))
   s <- summary(lm(as.formula(paste('prop ~ Disease_State +',paste(B_cov,collapse='+'))),data=d))$coefficients
   propdiff[[length(propdiff)+1]] <- data.table(organ='brain',version=v,cell=cc,beta_disease_SD=s['Disease_State',1],p=s['Disease_State',4],
     mean_prop_control=NA_real_,mean_prop_case=NA_real_)
   cors[[length(cors)+1]] <- data.table(organ='brain',version=v,cell=cc,r_with_red=cor(bd$blue,brain_props[[v]][,cc]),r_with_turquoise=cor(bd$turquoise,brain_props[[v]][,cc])) }
}
res <- rbindlist(res); propdiff <- rbindlist(propdiff); cors <- rbindlist(cors)
setnames(cors,c('r_with_red','r_with_turquoise'),c('r_with_AxisA_source_ME','r_with_AxisB_source_ME'))
dir.create(file.path(OUT,'results'),showWarnings=FALSE)
fwrite(base,file.path(OUT,'results/baseline_reproduced.tsv'),sep='\t')
fwrite(res,file.path(OUT,'results/adjusted_models.tsv'),sep='\t')
fwrite(propdiff,file.path(OUT,'results/proportion_case_control.tsv'),sep='\t')
fwrite(cors,file.path(OUT,'results/ME_proportion_correlations.tsv'),sep='\t')
fwrite(rbindlist(lung_markers),file.path(OUT,'results/lung_dtangle_markers.tsv'),sep='\t')
fwrite(rbindlist(brain_markers),file.path(OUT,'results/brain_BRETIGEA_markers.tsv'),sep='\t')
saveRDS(list(lung=lung_props,brain=brain_props),file.path(OUT,'results/estimated_proportions.rds'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'logs/R_sessionInfo.txt'))
print(res[,.(organ,module,role,version,adjusted_for,beta_base=round(beta_base,3),beta_adj=round(beta_adj,3),lo=round(lo_adj,3),hi=round(hi_adj,3),p_adj=signif(p_adj,3),pct=round(100*pct_change,1))])

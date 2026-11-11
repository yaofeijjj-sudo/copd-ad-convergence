suppressPackageStartupMessages(library(data.table))
root<-normalizePath(getwd());out<-file.path(root,'cell_localization/SYN01_Entrez_downstream_20261001');model<-file.path(out,'results/models/tables')
for(n in c('lung_broad_model_profile_weights.tsv','brain_broad_model_profile_weights.tsv','lung_matched_null_draws.tsv','brain_matched_null_draws.tsv','lung_myeloid_CR2.tsv','brain_subtype_CR2.tsv','lung_myeloid_fine_donor_profiles.tsv','supplementary24_brain_broad_donor_profiles.tsv'))if(file.exists(file.path(model,n)))stopifnot(file.copy(file.path(model,n),file.path(out,'source_data',n),overwrite=TRUE))
for(organ in c('lung','brain')){
 n<-if(organ=='lung')'lung_myeloid_donor_profiles.tsv' else 'brain_subtype_donor_profiles.tsv'
 if(!file.exists(file.path(model,n)))next
 d<-fread(file.path(model,n));cov<-d[n_cells>=20,.(n_donors=uniqueN(donor)),by=annotation];d<-d[n_cells>=20 & annotation %chin% cov[n_donors>=5]$annotation]
 d[,original_annotation:=annotation];fwrite(d,file.path(out,'source_data',if(organ=='lung')'figure_C_donor_scores.tsv' else 'figure_D_donor_scores.tsv'),sep='\t')
}
for(n in c('brain_subtype_display_map.tsv'))file.copy(file.path(root,'cell_localization/subtype_raincloud_color_20260917/source_data',n),file.path(out,'source_data',n),overwrite=TRUE)
file.copy(file.path(root,'cell_localization/independent_figures_20260917/source_data/brain_subtype_display_map_20260929.tsv'),file.path(out,'source_data/brain_subtype_display_map_20260929.tsv'),overwrite=TRUE)

# The 24 official subclasses are descriptive supplementary profiles, with the same mean score.
f<-file.path(out,'results/brain_work/broad_ucell_v1/profile_gene_rank_aggregate.rds')
if(file.exists(f)){a<-readRDS(f);p<-copy(a$profiles);p[,n_cells:=a$profile_n_cells];m<-length(a$axis_genes);p[,mean_ucell:=colSums(1500-a$mean_truncated_rank[a$axis_genes,,drop=FALSE])/(m*1500-m*(m+1)/2)];fwrite(p[n_cells>=20],file.path(model,'supplementary24_brain_broad_donor_profiles.tsv'),sep='\t');fwrite(p[n_cells>=20],file.path(out,'source_data/supplementary24_brain_broad_donor_profiles.tsv'),sep='\t')}

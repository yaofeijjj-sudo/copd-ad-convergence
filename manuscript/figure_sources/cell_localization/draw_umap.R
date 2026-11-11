suppressPackageStartupMessages({library(data.table);library(grid)})
root<-normalizePath(getwd());out<-file.path(root,'细胞定位/SYN01_Entrez_downstream_20261001')
args<-commandArgs(trailingOnly=TRUE);organs<-if(length(args))args else c('Lung','Brain')
for(organ in organs){
if(organ=='Brain'){
 src<-file.path(root,'细胞定位/brain_umap_reference_20260918/source_data/brain_umap_cells.tsv.gz')
 d<-fread(src);stopifnot(nrow(d)==1357145)
 z<-fread(file.path(out,'results/brain_work/broad_ucell_v1/cell_scores.tsv'),select=c('cell_id','axisA_full_measured_ucell'));stopifnot(nrow(z)==nrow(d),!anyDuplicated(z$cell_id));idx<-match(d$cell_id,z$cell_id);stopifnot(!anyNA(idx));d[,score:=z$axisA_full_measured_ucell[idx]]
 palette<-c('Immune'='#F04438','Astrocyte'='#F28C45','Endothelial'='#DDB34B','Excitatory neurons'='#52B7A5','Inhibitory neurons'='#4C9FC4','OPC'='#8571C4','Oligodendrocyte'='#D76E98','VLMC & Perivascular'='#8DA85D')
}else{
 src<-file.path(root,'细胞定位/lung_umap_reference_20260922/source_data/lung_umap_cells.tsv.gz')
 d<-fread(src);stopifnot(nrow(d)==211464)
 z<-fread(file.path(out,'results/lung/objects/axisA_cell_scores.tsv.gz'),select=c('cell_id','AxisA_UCell'));stopifnot(nrow(z)==nrow(d),!anyDuplicated(z$cell_id));idx<-match(d$cell_id,z$cell_id);stopifnot(!anyNA(idx));d[,`:=`(score=z$AxisA_UCell[idx],group=broad_celltype)]
 palette<-c('Myeloid'='#F04438','Epithelial'='#F28C45','Endothelial'='#DDB34B','T/NK'='#52B7A5','B/Plasma'='#8571C4','Stromal'='#D76E98','Other/Progenitor'='#AEB5BF')
}
stopifnot(all(d$group %in% names(palette)))
fwrite(data.table(group=names(palette),color=unname(palette)),file.path(out,'source_data',paste0(organ,'_palette.tsv')),sep='\t')
xr<-range(d$UMAP_1);yr<-range(d$UMAP_2);scale<-min(150/diff(xr),140/diff(yr));cx<-93;cy<-97
px<-cx+(d$UMAP_1-mean(xr))*scale;py<-cy+(d$UMAP_2-mean(yr))*scale
ramp<-colorRampPalette(c('#E2E5E9','#F7B0A1','#F04438'))(256)
stopifnot(all(is.finite(d$score)),all(d$score>=0),all(d$score<=1));fwrite(d[,.(cell_id,donor=if(organ=='Brain')donor else subject_id,group,UMAP_1,UMAP_2,AxisA_Entrez_UCell=score)],file.path(out,'source_data',paste0(organ,'_UMAP_cells.tsv.gz')),sep='\t');scoremax<-max(d$score);bins<-floor(d$score/scoremax*255)+1L
preview<-file.path(root,'_workspace/SYN01_downstream_update_20261001');dir.create(preview,recursive=TRUE,showWarnings=FALSE)
for(kind in c('FeaturePlot','Cell_types')){
 name<-paste0(organ,'_UMAP_',kind);labs<-character()
 quartz(type='pdf',file=file.path(out,'figures',paste0(name,'.pdf')),width=270/25.4,height=190/25.4,family='Arial',pointsize=10,bg='white');grid.newpage()
 T<-function(s,x,y,just='left'){grid.text(s,unit(x,'mm'),unit(y,'mm'),just=just,gp=gpar(fontfamily='Arial',fontsize=10,fontface='bold',col='#262626'));labs<<-c(labs,s)}
 T(paste0(organ,if(kind=='FeaturePlot')': Axis A UCell' else ': cell types'),18,180)
 if(kind=='FeaturePlot'){ord<-order(d$score,seq_len(nrow(d)));cols<-ramp[bins]}else{set.seed(18);ord<-sample.int(nrow(d));cols<-unname(palette[d$group])}
 # Every nucleus is retained; opaque vector circles, no rasterization. High scores drawn last.
 for(start in seq(1,length(ord),by=50000)){
  idx<-ord[start:min(start+49999,length(ord))]
  grid.circle(unit(px[idx],'mm'),unit(py[idx],'mm'),r=unit(.085,'mm'),gp=gpar(fill=cols[idx],col=NA))
 }
 grid.lines(unit(c(18,18,38),'mm'),unit(c(43,23,23),'mm'),gp=gpar(col='#262626',lwd=1.4))
 T('UMAP 1',18,15);T('UMAP 2',18,49)
 if(kind=='Cell_types'){
  for(i in seq_along(palette)){yy<-145-(i-1)*13;grid.circle(unit(179,'mm'),unit(yy,'mm'),r=unit(1.5,'mm'),gp=gpar(fill=palette[i],col=NA));T(names(palette)[i],184,yy)}
 }else{
  T('Axis A UCell',180,150)
  for(i in 1:256)grid.rect(unit(180,'mm'),unit(65+(i-1)*70/256,'mm'),unit(5,'mm'),unit(70/256+.01,'mm'),just=c('left','bottom'),gp=gpar(fill=ramp[i],col=NA))
  for(v in seq(0,.2,.05)[seq(0,.2,.05)<scoremax*.94])T(sprintf('%.2f',v),188,65+v/scoremax*70)
  T(sprintf('%.3f',scoremax),188,135)
 }
 dev.off();fwrite(data.table(label=labs),file.path(out,'qa',paste0(name,'_labels.tsv')),sep='\t')
 stopifnot(system2(Sys.which('pdftoppm'),c('-scale-to','1600','-png','-singlefile',shQuote(file.path(out,'figures',paste0(name,'.pdf'))),shQuote(file.path(preview,name))))==0)
 fwrite(data.table(organ=organ,kind=kind,n_cells=nrow(d),width_mm=270,height_mm=190,x_center_mm=cx,y_center_mm=cy,coordinate_scale=scale,radius_mm=.085,score_max=scoremax,coordinates_unchanged=TRUE),file.path(out,'qa',paste0(name,'_geometry.tsv')),sep='\t');message('Completed ',name)
}
writeLines(c(capture.output(sessionInfo()),'Full observed score range; no clipping or smoothing. Equal coordinate scales across both PDFs. Vector circle radius 0.085 mm. FeaturePlot ascending score order; categorical seed 18 shuffled order.'),file.path(out,'qa/drawing_specification.txt'))

}

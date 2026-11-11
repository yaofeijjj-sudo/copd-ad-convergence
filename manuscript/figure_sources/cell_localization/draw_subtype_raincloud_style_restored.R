suppressPackageStartupMessages({library(data.table);library(ggplot2);library(grid)})
root<-normalizePath(getwd());out<-file.path(root,'cell_localization/SYN01_Entrez_downstream_20261001')
draw_one<-function(organ,file,expected,n_groups,xmax,title,palette){
d<-fread(file.path(out,'source_data',file));d[,annotation:=original_annotation]
labelmap<-setNames(unique(d$annotation),unique(d$annotation))
if(organ=='Brain'){
 m<-fread(file.path(out,'source_data/brain_subtype_display_map.tsv'))
 m[,display_label:=gsub('\\n','\n',display_label,fixed=TRUE)]
 labelmap<-setNames(m$display_label,m$original_label)
}
if(organ=='Lung'){labelmap<-c('Alveolar macrophages'='Alveolar\nmacrophages','Monocyte-derived macrophages'='Monocyte-derived\nmacrophages','Perivascular interstitial macrophages'='Perivascular interstitial\nmacrophages','Classical monocytes'='Classical\nmonocytes','Conventional DCs'='Conventional\nDCs','Non-classical monocytes'='Non-classical\nmonocytes','Plasmacytoid DCs'='Plasmacytoid\nDCs','Mast cells'='Mast cells')}
stopifnot(nrow(d)==expected,uniqueN(d$annotation)==n_groups,anyDuplicated(d[,.(donor,annotation)])==0,capabilities('aqua'))
ord<-d[,.(med=median(mean_ucell)),by=annotation][order(-med)]$annotation
grDevices::quartz(type='pdf',file=tempfile(fileext='.pdf'),width=9,height=6,family='Arial',pointsize=14)
labelwidth<-max(vapply(unlist(strsplit(unname(labelmap[ord]),"\n",fixed=TRUE)),function(z)convertWidth(grobWidth(textGrob(z,gp=gpar(fontfamily='Arial',fontsize=14,fontface='bold'))),'mm',TRUE),numeric(1)));dev.off()
left<-max(48,labelwidth+9);right<-left+118;w<-right+10;bottom<-25;gap<-17;h<-bottom+18+n_groups*gap;row_y<-h-23-(seq_len(n_groups)-1)*gap
x<-function(v)left+v/xmax*(right-left)
font<-'Arial';fontsize<-14;ink<-'#262626';coral<-'#F04438';grey<-'#AEB5BF'
records<-list();pointdata<-list();densitydata<-list();boxdata<-list()
T<-function(s,x,y,just='left',bold=TRUE){
 lines<-strsplit(s,'\n',fixed=TRUE)[[1]]
 for(j in seq_along(lines)){
 yy<-y+((length(lines)+1)/2-j)*5.6
 grid.text(lines[j],unit(x,'mm'),unit(yy,'mm'),just=just,gp=gpar(fontfamily=font,fontsize=14,fontface=if(bold)'bold' else 'plain',col=ink))
 records[[length(records)+1]]<<-data.table(label=lines[j])
 }
}
L<-function(x1,y1,x2,y2,col=ink,lwd=1.4)grid.lines(unit(c(x1,x2),'mm'),unit(c(y1,y2),'mm'),gp=gpar(col=col,lwd=lwd))
grDevices::quartz(type='pdf',file=file.path(out,'figures',paste0(organ,'_subtype_raincloud.pdf')),width=w/25.4,height=h/25.4,family=font,pointsize=14,bg='white');grid.newpage()
T(title,left,h-8,bold=TRUE)
L(left,bottom,right,bottom)
for(v in seq(0,xmax,.05)){L(x(v),bottom,x(v),bottom-1.5);T(format(v,trim=TRUE),x(v),bottom-5,just='centre')}
T('Mean Axis A UCell',(left+right)/2,7,just='centre')
set.seed(17)
for(i in seq_along(ord)){
 a<-ord[i];z<-copy(d[annotation==a]);yy<-row_y[i];col<-unname(palette[a]);stopifnot(!is.na(col));fill<-adjustcolor(col,.45);pointcol<-adjustcolor(col,.7)
 # Official ggplot2 stat_ydensity computes the half-violin geometry.
 vd<-as.data.table(ggplot_build(ggplot(z,aes(x=factor(1),y=mean_ucell))+geom_violin(trim=TRUE,scale='width',bw='nrd0',adjust=1))$data[[1]])
 setorder(vd,y)
 # Official ggplot2 stat_boxplot: median, quartiles and 1.5-IQR whiskers.
 bd<-as.data.table(ggplot_build(ggplot(z,aes(x=factor(1),y=mean_ucell))+geom_boxplot(coef=1.5))$data[[1]])
 T(unname(labelmap[a]),left-4,yy,just='right')
 grid.polygon(unit(c(x(vd$y),rev(x(vd$y))),'mm'),unit(c(yy+1+vd$violinwidth*5,rep(yy+1,nrow(vd))),'mm'),gp=gpar(fill=fill,col=col,lwd=1.4))
 by<-yy-1.2
 L(x(bd$ymin),by,x(bd$ymax),by,col=col)
 L(x(bd$ymin),by-.55,x(bd$ymin),by+.55,col=col);L(x(bd$ymax),by-.55,x(bd$ymax),by+.55,col=col)
 grid.rect(unit(x(bd$lower),'mm'),unit(by-.7,'mm'),unit(x(bd$upper)-x(bd$lower),'mm'),unit(1.4,'mm'),just=c('left','bottom'),gp=gpar(fill='white',col=col,lwd=1.4))
 L(x(bd$middle),by-.7,x(bd$middle),by+.7)
 z[,point_y_mm:=yy-3.6+runif(.N,-.7,.7)];z[,point_x_mm:=x(mean_ucell)]
 grid.circle(unit(z$point_x_mm,'mm'),unit(z$point_y_mm,'mm'),r=unit(.42,'mm'),gp=gpar(fill=pointcol,col=NA))
 z[,category_order:=i];pointdata[[i]]<-z
 densitydata[[i]]<-vd[,.(annotation=a,value=y,density=density,scaled_width=violinwidth)]
 boxdata[[i]]<-bd[,.(annotation=a,lower_whisker=ymin,q1=lower,median=middle,q3=upper,upper_whisker=ymax)]
}
dev.off()
fwrite(rbindlist(records),file.path(out,'qa',paste0(organ,'_text_labels.tsv')),sep='\t')
fwrite(data.table(organ=organ,annotation=ord,display=unname(labelmap[ord]),color=unname(palette[ord])),file.path(out,'source_data',paste0(organ,'_palette.tsv')),sep='\t')
fwrite(data.table(width_mm=w,height_mm=h,row_gap_mm=gap,left_mm=left,right_mm=right,xmax=xmax),file.path(out,'qa',paste0(organ,'_geometry.tsv')),sep='\t')
fwrite(rbindlist(pointdata),file.path(out,'source_data',paste0(organ,'_plotted_points.tsv')),sep='\t')
fwrite(rbindlist(densitydata),file.path(out,'source_data',paste0(organ,'_density_geometry.tsv')),sep='\t')
fwrite(rbindlist(boxdata),file.path(out,'source_data',paste0(organ,'_boxplot_summary.tsv')),sep='\t')
writeLines(c(capture.output(sessionInfo()),capture.output(systemfonts::match_fonts('Arial'))),file.path(out,'qa',paste0(organ,'_R_session.txt')))

}
lung_palette<-c('Alveolar macrophages'='#F04438','Monocyte-derived macrophages'='#F28C45','Perivascular interstitial macrophages'='#DDB34B','Classical monocytes'='#52B7A5','Conventional DCs'='#4C9FC4','Non-classical monocytes'='#8571C4','Plasmacytoid DCs'='#D76E98','Mast cells'='#8DA85D')
brain_palette<-c('Micro-PVM_2'='#F04438','Micro-PVM_3-SEAAD'='#F28C45','Micro-PVM_2_3-SEAAD'='#DDB34B','Micro-PVM_1'='#52B7A5','Micro-PVM_2_1-SEAAD'='#8571C4','Lymphocyte'='#D76E98')
requested<-commandArgs(trailingOnly=TRUE)
if(length(requested)==0L || 'Lung' %in% requested) draw_one('Lung','figure_C_donor_scores.tsv',287,8,.3,'Lung: myeloid subtypes',lung_palette)
if(length(requested)==0L || 'Brain' %in% requested) draw_one('Brain','figure_D_donor_scores.tsv',244,6,.2,'Brain: immune subtypes',brain_palette)

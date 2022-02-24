.libPaths(c("/home/che82/R/x86_64-pc-linux-gnu-library/4.1","/usr/local/lib/R/site-library","/usr/lib/R/site-library","/usr/lib/R/library" ))

library(rdist)
library(SCORPIUS)
library(Rmagic)
library(pheatmap)
library(Linnorm)
#######on alignment##########
set.seed(1000)
library(reshape2)
Sys.setenv(RETICULATE_PYTHON = "/home/che82/miniconda3/bin/python")
library(reticulate)
#set pathon path#########hcf
path_to_python<-'/home/che82/miniconda3/bin/python'
use_python(path_to_python)
##########################
library(ManiNetCluster)
library(plyr)
library(RColorBrewer)
library(stringr)
library(limma)

library(ComplexHeatmap)
library(circlize)
#######load and processing data#########
source('func.r')



#####load prepare data from align######
load('bulk_start.RData')

rdata1 = rdata_human; meta1=meta_human
rdata2 = rdata_organoid; meta2=meta_organoid


###################
#regions=c('MFC','OFC','DFC','VFC','M1C','S1C','IPC','A1C','STC','ITC','V1C','HIP','AMY','STR','MD','CBC')
regions = unique(meta1$Region)
form.data1=rdata1[,meta1$Species=='Human' & meta1$Region %in% regions];form.meta1=meta1[meta1$Species=='Human' & meta1$Region %in% regions,]
form.data2=rdata2;form.meta2=meta2



form.meta1$time=form.meta1$Period
form.meta2$time=form.meta2$Days


#??????? find significance DE in any condition#######

deg_list1 = DE.list(form.data1,form.meta1)
deg_list2 = DE.list(form.data2,form.meta2)


#####MAIN: start align##########
######focus on genes of interest#########
##co-exp genes+HVG genes
sel.genes = intersect(intersect(deg_list1,deg_list2),unique(all.rec$gene))
## co-expression genes
#sel.genes = intersect(intersect(row.names(form.data1),row.names(form.data2)),unique(all.rec$gene))
##all genes
#sel.genes = intersect(row.names(rdata1),row.names(rdata2))
##highly variable genes only
#sel.genes = intersect(row.names(form.data1),row.names(form.data2))

#select expression based genes,log transform and reorder
sel.data1 = form.data1[row.names(form.data1) %in% sel.genes,]; sel.data1 = sel.data1[!duplicated(row.names(sel.data1)),]; sel.meta1=form.meta1  #on human data
sel.data2 = form.data2[row.names(form.data2) %in% sel.genes,]; sel.data2 = sel.data2[!duplicated(row.names(sel.data2)),]; sel.meta2=form.meta2  #on organoid data

#reorder
exp1 = sel.data1[order(row.names(sel.data1)),order(sel.meta1$time)];sel.meta1=sel.meta1[order(sel.meta1$time),]
exp2 = sel.data2[order(row.names(sel.data2)),order(sel.meta2$time)];sel.meta2=sel.meta2[order(sel.meta2$time),]

#######below needs to refine#############
ps.mat1 = t(exp1);ps.time1=sel.meta1$time
ps.mat2 = t(exp2);ps.time2=sel.meta2$time


#align and visual
algn_res = runMSMA_dtw(ps.mat1,ps.mat2)

df2 = algn_res[[3]]

###########some plots############
#df2$time = c(sel.meta1$time,as.numeric(as.factor(sel.meta2$time)))
df2$time = c(sel.meta1$time,sel.meta2$time)
########calculate pairwise distances between cells after MSMA#####
pair_dist = apply(df2[df2$data=='sample1',c(3:5)],1,function(x) {
	d = apply(df2[df2$data=='sample2',c(3:5)],1,function(y) eu.dist(x,y))
})
row.names(pair_dist)=row.names(ps.mat2)
colnames(pair_dist)=row.names(ps.mat1)
sim_dist = 1/(1+pair_dist)
#hmtp3-2nd align
cols1 = c(brewer.pal(9,'YlGnBu')[c(2:9)],brewer.pal(11,'BrBG')[c(8:11)]);names(cols1) = unique(ps.time1)
annot1=rowAnnotation(H_time=factor(ps.time1,levels=unique(ps.time1)),col=list(H_time=cols1))
cols2 = c(brewer.pal(9,'YlOrRd')[c(2:9)],brewer.pal(11,'BrBG')[rev(1:4)]);names(cols2) = unique(ps.time2)
annot2=columnAnnotation(O_time=factor(ps.time2,levels=unique(ps.time2)),col=list(O_time=cols2))

sim_mat = t(sim_dist)
Heatmap(sim_mat,name='human_vs_humanOrg',show_row_names=F,show_column_names=F,cluster_rows=F,cluster_columns=F,left_annotation=annot1,top_annotation=annot2,col=colorRamp2(c(0.8,max(sim_dist)),c('white','red')))



####plot 3D trajectors####
pdf('3d1.pdf')
time.cols1 = colorRampPalette(brewer.pal(n=9,'Greens'))(12)
#time.cols1 = c(brewer.pal(9,'YlGnBu')[c(2:9)],brewer.pal(11,'BrBG')[c(8:11)])
res = data.frame(df2[df2$data=='sample1',])
res0 = data.frame(df2)
library(plot3D) 
s3d<-scatter3D(x=res[,3],y=res[,4],z=res[,5],colvar=as.numeric(mapvalues(res$time,names(table(res$time)),c(2:13))),col = c(time.cols1),pch=c(16,17)[as.numeric(as.factor(res$data))],colkey=F,theta = 300, phi = 30,cex=2,
xlim=c(min(res0$Val0),max(res0$Val0)),
ylim=c(min(res0$Val1),max(res0$Val1)),
zlim=c(min(res0$Val2),max(res0$Val2)))
legend("top", legend = levels(as.factor(res$data)), pch = c(16, 17),inset = -0.1, xpd = TRUE, horiz = TRUE)
legend("right", legend = levels(as.factor(res$time)), col = c(time.cols1),pch=16,inset =-0.05, xpd = TRUE, horiz = F,cex=1.2)
dev.off()

#time.cols2 = colorRampPalette(brewer.pal(n=12,'Set2'))(12)
#time.cols2 = c(brewer.pal(9,'YlOrRd')[c(2:9)],brewer.pal(11,'BrBG')[rev(1:4)])
#c25 <- c( "dodgerblue2", "#E31A1C",  "green4", "#6A3D9A",  "#FF7F00",  "black", "gold1", "skyblue2", "#FB9A99",  "palegreen2", "#CAB2D6",  "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown")
#time.cols2 = c25[1:12]

#time.cols2 = colorRampPalette(brewer.pal(9, "YlOrBr"))(12) #for human
time.cols2 = colorRampPalette(brewer.pal(9, "Blues"))(12) #for organoid

res = data.frame(df2[df2$data=='sample2',])
#res = data.frame(df2)
res0 = data.frame(df2)
library(plot3D) 
pdf('3D2.pdf',width=10)
s3d<-scatter3D(x=res[,3],y=res[,4],z=res[,5],pch=24,col='gray',bg=time.cols2[as.numeric(mapvalues(res$time,names(table(res$time)),c(1:12)))],colkey=F,theta = 300, phi = 30,cex=2,lwd=1,
xlim=c(min(res0$Val0),max(res0$Val0)),
ylim=c(min(res0$Val1),max(res0$Val1)),
zlim=c(min(res0$Val2),max(res0$Val2)))
legend("top", legend = levels(as.factor(res$data)), pch = c(16, 17),inset = -0.1, xpd = TRUE, horiz = TRUE)
legend("right", legend = levels(as.factor(res$time)), col = c(time.cols2),pch=16,inset = -0.05, xpd = TRUE, horiz = F,cex=1.2)
dev.off()


#######corrplot on timepoint wise averaged similarity
library(corrplot)
sim_avg = matrix(0,nrow=length(unique(ps.time1)),ncol=length(unique(ps.time2)))

i = 0
for(t1 in unique(ps.time1)) {
	i = i+1
	j = 0
	for (t2 in unique(ps.time2)) {
		j = j+1
		sim_tmp = sim_mat[ps.time1==t1,ps.time2==t2]
		avg_tmp = mean(sim_tmp)
		sim_avg[i,j] = avg_tmp
	}
}
sim_avg = t(sim_avg)
colnames(sim_avg) = unique(ps.time1)
row.names(sim_avg) = unique(ps.time2)

library(RColorBrewer)
pdf('corr.plot.bulk.human.vs.org1.pdf')
corrplot(sim_avg,col=rev(colorRampPalette(brewer.pal(9,'PRGn'))(200)),cl.lim=c(0.7,1),is.corr=F,tl.col="black")
dev.off()






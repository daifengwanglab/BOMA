library(Seurat)
library(rdist)


#functions

DE.list<-function(data,meta) {
        genes = c()
        for (t in unique(meta$time)) {
                condition=as.factor(meta$time==t)
                design = model.matrix(~0+condition)  # Create a design matrix containing time information             
                fit=lmFit(data,design)  # Linear model fitting
                fit = eBayes(fit)  # Bayesian test  
                top_table = topTable(fit,number=10000)  # Get the differentially expressed gene name
                genes = c(genes,row.names(top_table[top_table$adj.P.Val<0.05,]))  # combine into a list
                 }
        genes = unique(genes)
        return(genes)
}

get_form_data<-function(rdata,meta,type='NORM',hvg=6000,add=NULL) {#prefiltering and formalize the expression data,type='NORM' or 'COUNTS'
	form.data = c()
	if (type=='COUNTS') {
		obj <- CreateSeuratObject(counts = rdata, min.cells = 30, min.features = 200)
		obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
		obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
		obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = hvg)
		form.data <- as.matrix(RelativeCounts(data = obj[['RNA']]@data,scale.factor=1e6))
	} else if (type=='NORM') {
		form.data <- rdata
	}
	if (is.null(add)) {
		form.data = form.data[row.names(form.data) %in% VariableFeatures(obj),]
		form.meta = meta[row.names(meta) %in% colnames(form.data),]
	} else {
		form.data = form.data[row.names(form.data) %in% unique(unlist(c(VariableFeatures(obj),add))),]
		form.meta = meta[row.names(meta) %in% colnames(form.data),]
	}
	return(list(form.data,form.meta))
}

eu.dist<-function(p1,p2) { #given 2 vectors of positions, calculate  distance[e.g.eucleadian]
	return(sqrt(sum((p1-p2)^2)))
}

inn.dist<-function(p1,p2) { #normalized inner product
	p1=p1/norm(p1,type="2")
	p2=p2/norm(p2,type="2")
	return(p1 %*% p2)
}


get_pseudo_cells<-function(gex,keep.perc=.2,tag=null,pcn=20) {#form/average cells into pseudo cells---original sc exps too stochastic,input is the log-transformed matrix; keep.perc: percentage of pseudo cells VS original cells(%)
	n = dim(gex)[1]
	pca.res = prcomp(gex)
	pcs = gex %*% pca.res$rotation
	tpc = pcs[,1:pcn] #top PCs only 
	d = rdist(tpc)
	psm = cutree(hclust(d),k=min(floor(n*keep.perc),n))  #membership of cell clusters
	
	#average expression of each cell cluster
	avg.gex = matrix(0,nrow=length(unique(psm)),ncol=dim(gex)[2])
	for (id in unique(psm)) {
		tmp = gex[psm==id,]
		if (!is.null(dim(tmp))) {
			avg.gex[id,] = colMeans(tmp)
		} else {
			avg.gex[id,] = tmp #need or not? throw singletons?
		}
	}
	tmp.tag = rep(tag,dim(avg.gex)[1])
	ps.names = paste(c(1:dim(avg.gex)[1]),tmp.tag,sep='-')
	row.names(avg.gex) = ps.names
	colnames(avg.gex) = colnames(gex)
	
	return(list(avg.gex,tmp.tag))
}


get_pseudo_cells2<-function(gex,keep.perc=.2,tag=null,pcn=20,ctp=null,ctp.collapse.thr=.8) {
#form/average cells into pseudo cells---original sc exps too stochastic,input is the log-transformed matrix; keep.perc: percentage of pseudo cells VS original cells(%)
#if cell types list provided, also collapse on the cell type of pseudocell
	n = dim(gex)[1]
	pca.res = prcomp(gex)
	ctp = as.character(ctp);ctp[is.na(ctp)]='NA'
	pcs = gex %*% pca.res$rotation
	tpc = pcs[,1:pcn] #top PCs only 
	d = rdist(tpc)
	psm = cutree(hclust(d),k=min(floor(n*keep.perc),n))  #membership of cell clusters
	
	#average expression of each cell cluster
	avg.gex = matrix(0,nrow=length(unique(psm)),ncol=dim(gex)[2])
	ps.ctp = rep('NA',length(unique(psm)))
	for (id in unique(psm)) {
		tmp = gex[psm==id,]
		tmp.ctp = ctp[psm==id]
		if (!is.null(dim(tmp))) {
			avg.gex[id,] = colMeans(tmp)
			sorts = sort(table(tmp.ctp)/sum(table(tmp.ctp)),decreasing=T)
			if (sorts[1] >= ctp.collapse.thr) {
				ps.ctp[id] = names(sorts)[1]
			}
		} else {#need or not? throw singletons?
			avg.gex[id,] = tmp 
			ps.ctp[id] = tmp.ctp
		}
	}
	tmp.tag = rep(tag,dim(avg.gex)[1])
	ps.names = paste(c(1:dim(avg.gex)[1]),tmp.tag,sep='-')
	row.names(avg.gex) = ps.names
	colnames(avg.gex) = colnames(gex)
	tmp.tag = data.frame('tag'=rep(tag,dim(avg.gex)[1]),'ctp'=ps.ctp)
	
	return(list(avg.gex,tmp.tag))
}


get_pseudo_cells3<-function(gex,keep.perc=.2,tag=null,pcn=20,ctp=null,ctp.collapse.thr=.8) {
#form/average cells into pseudo cells---original sc exps too stochastic,input is the log-transformed matrix; keep.perc: percentage of pseudo cells VS original cells(%)
#if cell types list provided, also collapse on the cell type of pseudocell
#record which cell belongs to which pseudocell
	n = dim(gex)[1]
	pca.res = prcomp(gex)
	ctp = as.character(ctp);ctp[is.na(ctp)]='NA'
	pcs = gex %*% pca.res$rotation
	tpc = pcs[,1:pcn] #top PCs only 
	d = rdist(tpc)
	psm = cutree(hclust(d),k=min(floor(n*keep.perc),n))  #membership of cell clusters
	
	#average expression of each cell cluster
	avg.gex = matrix(0,nrow=length(unique(psm)),ncol=dim(gex)[2])
	ps.ctp = rep('NA',length(unique(psm)))
	for (id in unique(psm)) {
		tmp = gex[psm==id,]
		tmp.ctp = ctp[psm==id]
		if (!is.null(dim(tmp))) {
			avg.gex[id,] = colMeans(tmp)
			sorts = sort(table(tmp.ctp)/sum(table(tmp.ctp)),decreasing=T)
			if (sorts[1] >= ctp.collapse.thr) {
				ps.ctp[id] = names(sorts)[1]
			}
		} else {#need or not? throw singletons?
			avg.gex[id,] = tmp 
			ps.ctp[id] = tmp.ctp
		}
	}
	tmp.tag = rep(tag,dim(avg.gex)[1])
	ps.names = paste(c(1:dim(avg.gex)[1]),tmp.tag,sep='-')
	row.names(avg.gex) = ps.names
	colnames(avg.gex) = colnames(gex)
	tmp.tag = data.frame('tag'=rep(tag,dim(avg.gex)[1]),'ctp'=ps.ctp)
	
	tmp.belong = paste(psm,tag,sep='-')
	
	return(list(avg.gex,tmp.tag,tmp.belong,gex))
}


runMSMA_cor<-function(mat1,mat2,k=5) {
#mat hasDE.list<-function(data,meta) {#find DE genes at any time and combine into a list
        genes = c()
        for (t in unique(meta$time)) {
                condition=as.factor(meta$time==t)
                design = model.matrix(~0+condition)
                fit=lmFit(data,design)
                fit = eBayes(fit)
                top_table = topTable(fit,number=10000)
                genes = c(genes,row.names(top_table[top_table$adj.P.Val<0.05,]))
        }
        genes = unique(genes)
        return(genes)
}
runMSMA_mw<-function(mat1,mat2,k=5){

#mat has genes on the column !!!!!!!!!!!!!!!!!
	sim_mat = c() #matrix to store the similarity results after 1st step
	knn_ini = c() #matrix to store the correspondence matrix after 1st step align
	df1 = NULL     #for validation

	sim_mat = cor(t(mat1),t(mat2))
	#knn
	knn_ini=t(matrix(0,nrow=dim(sim_mat)[2],ncol=dim(sim_mat)[1]))
	cm = apply(sim_mat,2,function(x) order(x,decreasing=T)[1:k])
	for (i in c(1:dim(cm)[2])) {
		knn_ini[cbind(cm[,i]),i] = 1
	}
	
	#after get the knn_ini, use it as correspondence matrix for 2nd align
	XY_corr=Correspondence(matrix=knn_ini)
	df2=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',corr=XY_corr,d=3L,method='nonlinear manifold aln',k_NN=3L,k_medoids=6L) #NMA
	
	return(list(sim_mat,knn_ini,df2,df1))

}


runMSMA_mw<-function(mat1,mat2,k=5) {
#mat has genes on the column !!!!!!!!!!!!!!!!!
	sim_mat = c() #matrix to store the similarity results after 1st step
	knn_ini = c() #matrix to store the correspondence matrix after 1st step align
	df1 = NULL     #for validation

	df=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',d=10L,method='manifold warping',k_NN=3L,k_medoids=3L)
	pair_dist = apply(df[df$data=='sample1',c(3:12)],1,function(x) {
			d = apply(df[df$data=='sample2',c(3:12)],1,function(y) eu.dist(x,y))
	})
	row.names(pair_dist)=row.names(mat2)
	colnames(pair_dist)=row.names(mat1)
	range01 <- function(x){(x-min(x))/(max(x)-min(x))}
	pair_dist = range01(pair_dist)
	sim_mat = 1/(1+pair_dist)
	df1 = df  #range01 OR NOT?
	#knn
	k1=k;k2=k
	knn_mat1 = matrix(0,nrow=dim(pair_dist)[1],ncol=dim(pair_dist)[2]) #for row
	knn_mat2 = matrix(0,nrow=dim(pair_dist)[1],ncol=dim(pair_dist)[2]) #for col
	for(i in 1:dim(knn_mat1)[1]) {#for row
		knn_mat1[i,order(pair_dist[i,])[1:k1]]=1
	}
	for(i in 1:dim(knn_mat2)[2]) {#for col
			knn_mat2[order(pair_dist[,i])[1:k2],i]=1
	}
	knn_ini = t(1*(knn_mat1 & knn_mat2))
	
	#after get the knn_ini, use it as correspondence matrix for 2nd align
	XY_corr=Correspondence(matrix=knn_ini)
	df2=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',corr=XY_corr,d=3L,method='nonlinear manifold aln',k_NN=3L,k_medoids=6L) #NMA
	
	return(list(sim_mat,knn_ini,df2,df1))

}


runMSMA_dtw<-function(mat1,mat2) { #method can be 'cor','mw' or 'dtw'
#mat has genes on the column !!!!!!!!!!!!!!!!!
#1st step has 3 options:1)correlation;2)Manifod warping;3)Dynamic Time Warping
	sim_mat = c() #matrix to store the similarity results after 1st step
	knn_ini = c() #matrix to store the correspondence matrix after 1st step align
	df1 = NULL     #for validation

	library(dtw)
	cor_mat = cor(t(mat1),t(mat2))
	dist_mat = 1/(1+cor_mat)
	dtw_res=dtw(dist_mat,keep=TRUE,step=asymmetric,open.end=T,open.begin=T)  #open beginning and ending for flexibility
	knn_ini = matrix(0,nrow=nrow(dist_mat),ncol=ncol(dist_mat))
	knn_ini[cbind(dtw_res$index1,dtw_res$index2)]=1
	sim_mat=NULL #position holder
	
	#after get the knn_ini, use it as correspondence matrix for 2nd align
	XY_corr=Correspondence(matrix=knn_ini)
	df2=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',corr=XY_corr,d=3L,method='nonlinear manifold aln',k_NN=3L,k_medoids=6L) #NMA
	
	return(list(sim_mat,knn_ini,df2,df1))

}



runMSMA_diag<-function(mat1,mat2) { #method can be 'cor','mw' or 'dtw'
#mat has genes on the column !!!!!!!!!!!!!!!!!
#1st step has 3 options:1)correlation;2)Manifod warping;3)Dynamic Time Warping
	sim_mat = c() #matrix to store the similarity results after 1st step
	knn_ini = c() #matrix to store the correspondence matrix after 1st step align
	df1 = NULL     #for validation

	if (nrow(mat1) != nrow(mat2)) {
		stop('for diag case---the numerb of samples need to be same')
	}else {
		knn_ini = diag(nrow(mat1))
	}
	
	#after get the knn_ini, use it as correspondence matrix for 2nd align
	XY_corr=Correspondence(matrix=knn_ini)
	df2=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',corr=XY_corr,d=3L,method='nonlinear manifold aln',k_NN=3L,k_medoids=6L) #NMA
	
	return(list(sim_mat,knn_ini,df2,df1))

}








runMSMA<-function(mat1,mat2,method='cor',k=5) { #method can be 'cor','mw' or 'dtw'
#mat has genes on the column !!!!!!!!!!!!!!!!!
#1st step has 3 options:1)correlation;2)Manifod warping;3)Dynamic Time Warping
	sim_mat = c() #matrix to store the similarity results after 1st step
	knn_ini = c() #matrix to store the correspondence matrix after 1st step align
	df1 = NULL     #for validation
	if (method=='cor') {
		sim_mat = cor(t(mat1),t(mat2))
		
		#knn
		knn_ini=t(matrix(0,nrow=dim(sim_mat)[2],ncol=dim(sim_mat)[1]))
		cm = apply(sim_mat,2,function(x) order(x,decreasing=T)[1:k])
		for (i in c(1:dim(cm)[2])) {
			knn_ini[cbind(cm[,i]),i] = 1
		}
	}else if (method=='mw') {
		df=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',d=10L,method='manifold warping',k_NN=3L,k_medoids=3L)
		pair_dist = apply(df[df$data=='sample1',c(3:12)],1,function(x) {
				d = apply(df[df$data=='sample2',c(3:12)],1,function(y) eu.dist(x,y))
		})
		row.names(pair_dist)=row.names(mat2)
		colnames(pair_dist)=row.names(mat1)
		range01 <- function(x){(x-min(x))/(max(x)-min(x))}
		pair_dist = range01(pair_dist)
		sim_mat = 1/(1+pair_dist)
		df1 = df  #range01 OR NOT?
		#knn
		k1=k;k2=k
		knn_mat1 = matrix(0,nrow=dim(pair_dist)[1],ncol=dim(pair_dist)[2]) #for row
		knn_mat2 = matrix(0,nrow=dim(pair_dist)[1],ncol=dim(pair_dist)[2]) #for col
		for(i in 1:dim(knn_mat1)[1]) {#for row
			knn_mat1[i,order(pair_dist[i,])[1:k1]]=1
		}
		for(i in 1:dim(knn_mat2)[2]) {#for col
				knn_mat2[order(pair_dist[,i])[1:k2],i]=1
		}
		knn_ini = t(1*(knn_mat1 & knn_mat2))

	}else if (method=='dtw') {
		library(dtw)
		cor_mat = cor(t(mat1),t(mat2))
		dist_mat = 1/(1+cor_mat)
		dtw_res=dtw(dist_mat,keep=TRUE,step=asymmetric,open.end=T,open.begin=T)  #open beginning and ending for flexibility
		knn_ini = matrix(0,nrow=nrow(dist_mat),ncol=ncol(dist_mat))
		knn_ini[cbind(dtw_res$index1,dtw_res$index2)]=1
		sim_mat=NULL #position holder
	}else if (method=='diag') {#use diagnal matrix as correspondence matrix---for 1-1 match
		if (nrow(mat1) != nrow(mat2)) {
			stop('for diag case---the numerb of samples need to be same')
		}else {
			knn_ini = diag(nrow(mat1))
		}
	}else {
		stop('not recognized method')
	}
	
	#after get the knn_ini, use it as correspondence matrix for 2nd align
	XY_corr=Correspondence(matrix=knn_ini)
	df2=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',corr=XY_corr,d=3L,method='nonlinear manifold aln',k_NN=3L,k_medoids=6L) #NMA
	
	return(list(sim_mat,knn_ini,df2,df1))

}









##########FUNCTION that identifies DEGs across two separate experiments#########











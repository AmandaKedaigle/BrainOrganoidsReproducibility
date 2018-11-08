#Code altered from @seanken

dir = '../round5_data'
outdir = '../round5_results/'

library(Seurat)
library(stringr)
library(Matrix)

#loads all 10X lanes from a given directory
dir10X<-function(dir="",outdir="",dat=NULL,lst=c(),makeSeurat=T,minGenes=500,regress=c("nUMI"),num=c())
{
	if(length(lst)==0) {	
		print(paste("ls ",dir,"/*/*/outs/filt*/* | grep : | sed 's/://g'",sep=""))
		lst=system(paste("ls ",dir,"/*/*/outs/filt*/* | grep : | sed 's/://g'",sep=""),intern=T)
	}	
	if(length(num)>0){
		lst=lst[num]
	}
	print(lst)
	if(is.null(dat)){
		print("Read in!")
		dat=Read10X(lst)
		print(head(dat))
		print("Fix colnames!")
		cols=colnames(dat)
		cols_new=c()
		for(col in cols){
			start=str_sub(col,1,1)
			cur=col
			if(start %in% c("A","T","G","C")){cur=paste("1_",cur,sep="")}
			cols_new<-c(cols_new,cur)
		}

		colnames(dat)=cols_new
		#print("Return!")
		saveRDS(dat, file = paste0(outdir,"dat.rds"))
	}
	print(paste("Dims: ",toString(dim(dat))))
	if(!makeSeurat){return(dat)}

	print("Make object!")
	seur<-CreateSeuratObject(dat,"Seurat",min.genes=minGenes,normalization.method="LogNormalize",scale.factor=1000000)


	#QC
	mito.genes <- grep(pattern = "^MT-", x = rownames(x = seur@data), value = TRUE)
	percent.mito <- Matrix::colSums(seur@raw.data[mito.genes, ]) / Matrix::colSums(seur@raw.data)
	seur <- AddMetaData(object = seur, metadata = percent.mito, col.name = "percent.mito")
	ribo.genes <- grep("^RP[S,L]",rownames(seur@raw.data), value = TRUE)
	percent.ribo <- Matrix::colSums(seur@raw.data[ribo.genes, ]) / Matrix::colSums(seur@raw.data)
	seur <- AddMetaData(object = seur, metadata = percent.ribo, col.name = "percent.ribo")
	pdf(paste0(outdir,'QC_Rplots.pdf'))
	par(mfrow = c(1, 3))
	GenePlot(object = seur, gene1 = "nUMI", gene2 = "percent.mito")
	GenePlot(object = seur, gene1 = "nUMI", gene2 = "percent.ribo")
	GenePlot(object = seur, gene1 = "nUMI", gene2 = "nGene")
	dev.off()

	print("Get variable genes!")
	par(mfrow = c(1,1))
	seur<-FindVariableGenes(seur,x.low.cutoff=1,do.plot =T)

	print("Regress out!")
	seur<-ScaleData(seur,genes.use=seur@var.genes,vars.to.regress=regress)

	print("Run PCA!")
	seur<-RunPCA(seur,pcs.compute=60)

	print("Save object!")
	saveRDS(seur, file = paste0(outdir,"initial_seur.rds"))

	pdf(paste0(outdir,'pca_Rplots.pdf'))
	VizPCA(object = seur, pcs.use = 1:2)
	PCAPlot(object = seur, dim.1 = 1, dim.2 = 2)
	PCHeatmap(object = seur, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
	PCElbowPlot(object = seur, 30)
	dev.off()
  return seur
}

seur = dir10X(dir=dir, outdir=outdir)

numPCs = 11 #change based on PCElbowPlot

library(Matrix)
library(reshape)
library(igraph)
library(RANN)

clust_Graph<-function(X,nn,getGraph=FALSE,type="louvain")
{

	nearest=nn2(X,X,k=nn+1, treetype = "bd", searchtype="priority")
	print("Found nearest neighbors")
	nearest$nn.idx = nearest$nn.idx[,-1]
	nearest$nn.dists = nearest$nn.dists[,-1] #Convert to a similarity score
	nearest$nn.sim = 1*(nearest$nn.dists >= 0 )
    
	edges = melt(t(nearest$nn.idx)); colnames(edges) = c("B", "A", "C"); edges = edges[,c("A","B","C")]
	edges$B = edges$C; edges$C=1

	edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))

	NN = nearest$nn.idx
	jaccard_dist = apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
           
	edges$C = jaccard_dist
	edges = subset(edges, C != 0)
	edges$C = edges$C/max(edges$C)
	colnames(edges) = c("First_Node", "Second_Node", "Weight")
	print("Make Sparse Adj!")
	adj_sparse<-sparseMatrix(i=edges$First_Node,j=edges$Second_Node,x=edges$Weight,symmetric=TRUE)
	print("Make graph!")
	g=graph_from_adjacency_matrix(adj_sparse,mode="undirected",weighted=TRUE)

	if(getGraph){return(g)}

	print("Cluster!")
	if(type=="louvain")
	{
	    clust = cluster_louvain(g)
	}
	else
	{
	    clust = cluster_infomap(g)
	}
	return(clust)
}

X=data.frame(seur@dr$pca@cell.embeddings)[1:numPCs]
clust<-clust_Graph(X,50,type="louvain")
seur@meta.data[paste("Clusts_nn","louvain",toString(50),"PC",toString(numPCs),sep="_")]=clust$membership;
seur<-SetAllIdent(seur,paste("Clusts_nn","louvain",toString(50),"PC",toString(numPCs),sep="_"))

print("Make tSNE!")
seur <- RunTSNE(object = seur, dims.use = 1:numPCs)
TSNEPlot(object = seur)
saveRDS(seur, file = paste0(outdir,'clusteredSeur.rds'))

seur.markers <- FindAllMarkers(object = seur, only.pos = T, test.use = "MAST", logfc.threshold = 0.25)
write.table(seur.markers, file="markerLists.txt")

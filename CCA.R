library(Seurat)
library(Matrix)
library(reshape)
library(igraph)
library(RANN)

#Provide seurat objects from individual datasets as arguments
args = commandArgs(trailingOnl=T)
outdir = args[1]
data1 = args[2]
data2 = args[3]
data3 = args[4]
dir.create(outdir)
setwd(outdir)

d1 = readRDS(data1)
d1@meta.data$dataset <- data1
d2 = readRDS(data2)
d2@meta.data$dataset <- data2
d3 = readRDS(data3)
d3@meta.data$dataset <- data3
print("Loaded data")

# Determine genes to use for CCA, must be highly variable in at least 2 datasets
ob.list <- list(d1,d2,d3)
genes.use <- c()
for (i in 1:length(ob.list)) {
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}

# Run multi-set CCA
combined <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 20)

# Run rare non-overlapping filtering
combined <- CalcVarExpRatio(object = combined, reduction.type = "pca",
                                       grouping.var = "dataset", dims.use = 1:10)
combined <- SubsetData(combined, subset.name = "var.ratio.pca",accept.low = 0.5)

#Align subspaces
numCCs = 20
combined <- AlignSubspace(combined, reduction.type = "cca", grouping.var = "dataset", 
                                 dims.align = 1:numCCs)

#tsne and clustering
combined <- RunTSNE(combined, reduction.use = "cca.aligned", dims.use = 1:numCCs, 
                           do.fast = T)

#ClustGraph function from readData file
X=data.frame(combined@dr$cca.aligned@cell.embeddings)[1:numCCs]
clust<-clust_Graph(X,50,type="louvain")
combined@meta.data[paste("Clusts_nn","louvain",toString(50),"CC",toString(numCCs),sep="_")]=clust$membership;
combined<-SetAllIdent(combined,paste("Clusts_nn","louvain",toString(50),"CC",toString(numCCs),sep="_"))

#Save results
combined@meta.data$org.num = sapply(strsplit(combined@cell.names, "_"), "[[",2)
combined@meta.data$org.ident = paste(combined@meta.data$dataset, combined@meta.data$org.num, sep=".")
saveRDS(combined, "cca_finished.rds")

#tsne plot each organoid
numOrgs = 9
combined <- SetAllIdent(combined, id="org.ident")
color.list=rep(c("#648FFF","#FFB000","#DC267F"),3)
p = list()
idents = c("1 PGP1.1","1 PGP1.2","1 HUES66", "2 PGP1.1","2 PGP1.2", "2 HUES66","3 PGP1.1","3 PGP1.2","3 HUES66")
for (num in 1:numOrgs){
  p[[num]] = TSNEPlot(combined, pt.size = 1, do.label = F, cells.use = WhichCells(combined, idents[num]),
                        no.axes=TRUE, no.legend=TRUE, colors.use = color.list[num])
  }
plot_grid(plotlist= p, ncol=3)

#tsne plot each line
combined <- SetAllIdent(combined,id="CellType")
combined@ident = factor(combined@ident, levels= c("Immature PNs","Immature CPNs","CPNs","Immature CFuPNs","CFuPNs", "IPCs/Immature PNs", "IPCs",
                                                  "Immature Interneurons","Ventral Precursors","Astroglia", "RG","oRG", 
                                                  "oRG/Astroglia","Cycling","Cajal-Retzius","Unknown"))
allcols = c("#332288","#88ccee", "#44aa99","#bbcc33","#117733","#999933","#999933","#ddcc77","#ee8866","#cc6677","#882255","#aa4499","#a51703","#c2a5cf","black",gray.colors(6)[6]) #From https://personal.sron.nl/~pault/
p = list()
cols = allcols[levels(combined@ident) %in% combined@ident[combined@meta.data$dataset=="PGP1.1"]]
p[[1]] = TSNEPlot(combined,no.axes=T,pt.size=1.5,no.legend=T, do.label=F, colors.use=cols, cells.use = combined@cell.names[combined@meta.data$dataset=="PGP1.1"])
cols = allcols[levels(combined@ident) %in% combined@ident[combined@meta.data$dataset=="PGP1.2"]]
p[[2]] = TSNEPlot(combined,no.axes=T,pt.size=1.5,no.legend=T, do.label=F, colors.use=cols, cells.use = combined@cell.names[combined@meta.data$dataset=="PGP1.2"])
cols = allcols[levels(combined@ident) %in% combined@ident[combined@meta.data$dataset=="HUES66"]]
p[[3]] = TSNEPlot(combined,no.axes=T,pt.size=1.5,no.legend=T, do.label=F, colors.use=cols, cells.use = combined@cell.names[combined@meta.data$dataset=="HUES66"])
plot_grid(plotlist= p, ncol=3)

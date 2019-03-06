library(Seurat)

numOrgs = 3
seur = readRDS('clusteredSeur.rds')
clustcolname = names(seur@meta.data)[[grep('^Clusts', names(seur@meta.data))]]

# Make the tsne plot
seur <- SetAllIdent(seur, id=clustcolname)
TSNEPlot(object = seur,do.label=T)

# Plot each organoid tsne separately
seur <- SetAllIdent(seur, id="orig.ident")
color.list = c("#ee8026","#dc050c","#72190e") #orange-brown
p = list()
p[[1]] = TSNEPlot(seur, pt.size = 0.3, do.label = F, colors.use = color.list, 
                  no.axes=TRUE, no.legend=TRUE)
for (num in 1:3){
  p[[num+1]] = TSNEPlot(seur, pt.size = 0.3, do.label = F, cells.use = WhichCells(seur, num), 
                        colors.use = color.list[num], 
                        no.axes=TRUE, no.legend=TRUE)
}
plot_grid(plotlist= p, ncol=4)

# Plot a tsne to look for background sources of variation in the data that could be affecting clustering
plot.by = c("nUMI","nGene","percent.mito","percent.ribo")
FeaturePlot(seur, features.plot = plot.by, pt.size = 0.5, no.legend = F)

#Calculate DEGs per cluster. Use this list to assign a cell type to each cluster
seur <- SetAllIdent(seur, id=clustcolname)
seur.markers <- FindAllMarkers(object = seur, only.pos = T, test.use = "MAST", logfc.threshold = 0.25)
write.table(seur.markers, file="markerLists.txt")

# tsne colored by cell type
# by now, should have ID'd cell types from DEG lists, and added a "CellType" column to the seur metadata
seur <- SetAllIdent(seur,id="CellType")
seur@ident = factor(seur@ident, levels= c("Immature PNs","Immature CPNs","CPNs","Immature CFuPNs","CFuPNs", "IPCs/Immature PNs", "IPCs",
                                          "Immature Interneurons","Ventral Precursors","Astroglia", "RG","oRG", 
                                          "Cycling","oRG/Astroglia","Unknown","Cajal-Retzius"))
cols = c("#332288","#88ccee", "#44aa99","#bbcc33","#117733", "#999933", "#999933", 
         "#ddcc77","#ee8866","#cc6677","#882255","#aa4499",
         "#c2a5cf","#a51703","darkgray","black") #From https://personal.sron.nl/~pault/
cols = cols[levels(seur@ident) %in% seur@ident]
TSNEPlot(seur,no.axes=T,pt.size=0.1,colors.use=cols,no.legend=T)

#Violin plots for some marker genes per celltype
g.lists = list(c("EOMES","PPP1R17","TMEM158"),
               c("HOPX","TNC","LGALS3"),
               c("SATB2","INHBA","FRMD4B"),
               c("BCL11B","CRYM","TLE4"),
               c("MKI67","TOP2A","BIRC5"),
               c("DLX1","DLX2","GAD2"))
types = c("IPC","oRG","CPN","CFuPN","Cycling","IN")
dir.create("VlnPlotsCellType")
for (t in 1:length(types)){
  png(paste0("VlnPlotsCellType/",types[t],'.vln.plots.png'),width = 1500, height=700, res=180)
  VlnPlot(seur, features.plot = g.lists[[t]], nCol=3, point.size.use=0, x.lab.rot = T,size.x.use = 10, 
                y.log=F, ident.include = c("Immature PNs","CPNs","CFuPNs", "Immature Interneurons", "IPCs","RG","oRG","Cycling"))
  dev.off()
}

# Percentage of cells in each cell type per organoid
seur <- SetAllIdent(seur, id="orig.ident")
library(reshape2)
counts = as.matrix(table(seur@meta.data[["CellType"]],seur@ident),2)
counts = melt(counts,varnames = c('CellType','Org'))
counts$CellType = factor(counts$CellType, levels= c("Immature PNs","Immature CPNs","CPNs","Immature CFuPNs","CFuPNs", "IPCs/Immature PNs", "IPCs",
                                               "Immature Interneurons","Ventral Precursors","Astroglia", "RG","oRG", 
                                               "oRG/Astroglia","Cycling","Unknown"))
cols = c("#332288","#88ccee", "#44aa99","#bbcc33","#117733", "#999933", "#999933", "#ddcc77","#ee8866","#cc6677","#882255","#aa4499","#a51703","#c2a5cf",gray.colors(6)[6]) #From https://personal.sron.nl/~pault/
cols = cols[levels(counts$CellType) %in% counts$CellType]
counts$Org = paste0("Org",as.numeric(counts$Org))
ggplot(counts, aes(fill=CellType, y=value, x=Org)) + 
  geom_bar(stat="identity", position="fill", show.legend=F)+
  labs(y="Percentage of Cells", x = "Sample") +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=16))+
  scale_fill_manual("legend",values=cols)

#Intraclass correlation between proportions of cell types in each organoid
res = data.frame(matrix(nrow=0, ncol=11))
Types = as.matrix(table(seur@meta.data[['CellType']],seur@meta.data$orig.ident),2)
for (ct in c("Immature PNs","Immature CPNs","CPNs","Immature CFuPNs","CFuPNs", "IPCs/Immature PNs", "IPCs",
             "Immature Interneurons","Ventral Precursors","Astroglia", "RG","oRG", 
             "Cycling","oRG/Astroglia","Unknown","Cajal-Retzius")) {
  if (ct %in% rownames(Types)){
    pNum = Types[ct,]/colSums(Types)
  } else {
    pNum=rep(0,ncol(Types))
  }
  res[ct,] = pNum
}
res["RG",] = res["oRG",] + res["oRG/Astroglia",] + res["RG",]
res = res[rownames(res)!="oRG/Astroglia" & rownames(res)!="oRG",]
res["CPNs",] = res["CPNs",] + res["Immature CPNs",]
res = res[rownames(res) != "Immature CPNs",]
resres["CFuPNs",] = res["CFuPNs",] + res["Immature CFuPNs",]
res = res[rownames(res) != "Immature CFuPNs",]
library(irr)
icc = icc(res, model="twoway",type="agreement", unit="single")
icc

# Pseudotime - in the end, ran this on CCA'd data from all batches together for each time point
library(monocle)
options(warn=-1)
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
                                           (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
pd <- new("AnnotatedDataFrame", data = seur@meta.data)
fd = data.frame("gene_short_name" = rownames(seur@raw.data))
rownames(fd) = rownames(seur@raw.data)
fd <- new("AnnotatedDataFrame", data = fd)
cds <- newCellDataSet(as.matrix(seur@raw.data),
                       phenoData = pd, featureData = fd,
                       expressionFamily=negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cell_type_color <- c("Immature PNs" = "#332288",
                     "Immature CPNs" = "#88ccee",
                     "CPNs" = "#44aa99",
                     "Immature CFuPNs" = "#bbcc33",
                     "CFuPNs" = "#117733",
                     "IPCs" = "#999933",
                     "IPCs/Immature PNs" = "#999933",
                     "Immature Interneurons" = "#ddcc77",
                     "Ventral Precursors" = "#ee8866",
                     "Astroglia" = "#cc6677",
                     'RG' = "#882255",
                     "oRG" = "#aa4499",
                     "Cycling" = "#c2a5cf",
                     "oRG/Astroglia" = "#a51703",
                     "Unknown" = "darkgray")
cds <- preprocessCDS(cds, num_dim = 20)
cds <- reduceDimension(cds, reduction_method = 'tSNE',verbose = T)
cds <- partitionCells(cds)
cds <- learnGraph(cds,  RGE_method = 'SimplePPT')
plot_cell_trajectory(cds,color_by = "CellType",show_tree = T) +
  scale_color_manual(values = cell_type_color) + theme(legend.position="none",axis.title = element_blank(),axis.text = element_blank(),
                                                       axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank())
MPP_node_ids = get_correct_root_state(cds,cell_phenotype ='CellType', "Cycling")
cds <- orderCells(cds, root_pr_nodes = MPP_node_ids)
plot_cell_trajectory(cds,show_tree = F) + theme(legend.position="none",axis.title = element_blank(),axis.text = element_blank(),
                                                axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank())


#Mutual information between clusters and organoids
#Downsample and recluster object
seur = SubsetData(combined,max.cells.per.ident=659)
seur<-FindVariableGenes(seur,x.low.cutoff=1,do.plot =T)
seur<-RunPCA(seur,pcs.compute=60)
PCElbowPlot(object = seur, 30)
numPCs = 20
X=data.frame(seur@dr$pca@cell.embeddings)[1:numPCs]
clust<-clust_Graph(X,50,type="louvain")
seur@meta.data$NewClusts=clust$membership;
#Create bar charts for cell clusters per organoid
library(reshape)
library(viridisLite)
counts = as.matrix(table(seur@meta.data[["NewClusts"]],seur@ident),2)
counts = melt(counts,varnames = c('Cluster','Sample'))
counts$Cluster = factor(as.character(counts$Cluster), levels=as.character(1:100))
counts$Sample = as.character(counts$Sample)
cols = viridis(max(seur@meta.data$NewClusts))
ggplot(counts, aes(fill=Cluster, y=value, x=Sample)) + 
  geom_bar(stat="identity", position="fill", show.legend=F)+
  scale_fill_manual(values = cols)+
  labs(y="Percentage of Cells", x = "Sample") +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=16))
#Calculate mutual information null model
save_dists = data.frame()
counts = seur@meta.data[,c("NewClusts","orig.ident")]
real_mi = dmi(counts)$mi[1,2]
randomized_mi = list()
for (i in 1:1000){
 randomClusts = sample(counts$NewClusts)
 counts$NewClusts = randomClusts
 randomized_mi[i] = as.numeric(dmi(counts)$mi[1,2])
}
randomized_mi = as.numeric(randomized_mi)
z = (real_mi - mean(randomized_mi)) / sd(randomized_mi)
#print MI score and z-score
real_mi
z

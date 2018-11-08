library(Seurat)
library(Matrix)
library(reshape)
library(igraph)
library(RANN)

data1 = "3mon_results/seur.RData"
data2 = "6mon_results/seur.RData"
outdir = "compare_3v6/"
setwd(outdir)

load(data1)
wt = seur
wt@meta.data$cond <- "3mon"
load(data2)
mut = seur
rm(seur)
mut@meta.data$cond <- "6mon"

# Gene selection for input to CCA
g.1 <- head(rownames(wt@hvg.info), 1000)
g.2 <- head(rownames(mut@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(wt@scale.data))
genes.use <- intersect(genes.use, rownames(mut@scale.data))

# Perform CCA (make single object and store vectors that project each dataset into the maximally correlated subspaces)
combined <- RunCCA(wt, mut, genes.use = genes.use, num.cc=30,add.cell.id1="3mon",add.cell.id2="6mon")

MetageneBicorPlot(combined, grouping.var = "stim", dims.eval = 1:30, display.progress = FALSE)
numCCs = 20

# Align subspaces
combined <- AlignSubspace(combined, reduction.type = "cca", grouping.var = "cond", dims.align = 1:numCCs)

# TSNE plot
combined <- RunTSNE(combined, reduction.use = "cca.aligned", dims.use = 1:numCCs, 
                           do.fast = T)
TSNEPlot(combined, do.return = T, pt.size = 0.5, group.by = "cond")

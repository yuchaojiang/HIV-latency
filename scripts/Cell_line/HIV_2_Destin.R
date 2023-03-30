setwd("~/Dropbox/HIV")

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(DescTools)
library(destin)
# https://api.github.com/repos/urrutiag/destin/tarball/HEAD
# install.packages("~/Desktop/urrutiag-destin-da94b69/destin.tar.gz", repos = NULL, type = "source")

load('processed_data/hiv.rda')

# Construct rse object
atac=hiv@assays$ATAC@counts
dim(atac)
rownames(atac)

atac[atac > 1] = 1
colData=colnames(atac)
rowData=rownames(atac)
rowRanges=StringToGRanges(rowData)
rse = SummarizedExperiment::SummarizedExperiment(assays = list(counts = atac),
                                                 rowRanges = rowRanges, colData = colData)

colData(rse)
rowRanges(rse)
assay(rse)[1:10,1:10]

# peak annotation
model = "hg38" # choose from hg19, hg38, mm10
rse = annotateRSE(rse, model)

# quality control
rse = doQC(rse, regionSumCutoff = 5, cellSumCutoffSDs = 3)

# Estimate number of clusters
clusterEst = estimateNClusters(rse, nClustersRange = 2:20)
nClusters = clusterEst$nClustersList$logLikeElbow
plotNClusters(clusterEst)
nClusters = 4

# Perform destin clustering
nCores = 2
clusterResults = destinGrid (rse, nClusters = nClusters, nCores = nCores)
clusterResults$cluster
cluster.results=clusterResults$cluster$cluster

# Visualization by tSNE
PCs = clusterResults$PCs
tsne = Rtsne(as.matrix(PCs))
componentX = tsne$Y[, 1]
componentY = tsne$Y[, 2]

toplot=data.frame(tsne1=componentX, tsne2=componentY,
                  group=hiv$group[names(cluster.results)])
sp<-ggplot(toplot, aes(x=tsne1, y=tsne2, color=group)) + geom_point(size = 0.5)+
  guides(colour = guide_legend(override.aes = list(size=5)))
sp

# Visualization by UMAP
library(reticulate)
use_python('/Users/yuchaojiang/pythonProject/bin/python', required = TRUE)

k=2; distance = "euclidean"; n_neighbors = 10; min_dist = 0.1; rand.seed = 42
set.seed(rand.seed)
reticulate::py_set_seed(rand.seed)
UMAP <- reticulate::import("umap")
umapper <- UMAP$UMAP(n_components = as.integer(k), metric = distance, 
                     n_neighbors = as.integer(n_neighbors), min_dist = min_dist)
Rumap <- umapper$fit_transform
umap <- Rumap(as.matrix(PCs))
componentX = umap[, 1]
componentY = umap[, 2]

toplot=data.frame(destinUMAP1=componentX, destinUMAP2=componentY,
                  group=hiv$group[names(cluster.results)])
sp<-ggplot(toplot, aes(x=destinUMAP1, y=destinUMAP2, color=group)) + geom_point(size = 0.5)+
  guides(colour = guide_legend(override.aes = list(size=5)))
sp

ggsave(sp, file='output/supp_figure1_destin.pdf', width=5, height=4)


DimPlot(hiv, reduction = "umap.atac", group.by = "group", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("scATAC-seq group")

toplot=data.frame(signacUMAP1=hiv@reductions$umap.atac@cell.embeddings[,1],
                  signacUMAP2=hiv@reductions$umap.atac@cell.embeddings[,2],
                  group=hiv$group)
sp<-ggplot(toplot, aes(x=signacUMAP1, y=signacUMAP2, color=group)) + geom_point(size = 0.5)+
  guides(colour = guide_legend(override.aes = list(size=5)))
sp

ggsave(sp, file='output/supp_figure1_signac.pdf', width=5, height=4)



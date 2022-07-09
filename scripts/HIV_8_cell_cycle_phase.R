setwd("~/Dropbox/HIV")
setwd("C:/Users/yuchaoj/Dropbox/HIV")

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

library(glmnet)
library(pheatmap)
library(fields)
library(qvalue)
library(gplots)
library(patchwork)
library(olsrr)
library(pdftools)
library(ape)
library(dendextend)

load('processed_data/hiv_predict_DE.rda')
hiv$pRNA=pRNA
hiv$pRNA.sqrt=sqrt(pRNA)

DimPlot(hiv, group.by='seurat_clusters', reduction='umap.rna', label = TRUE)
VlnPlot(hiv, features = 'pRNA.sqrt', group.by='seurat_clusters')
FeaturePlot(hiv, features = 'pRNA.sqrt', reduction = "umap.rna")

scores=readRDS(file='github/derived_data/cellcycle.aggr.rds')
all(colnames(hiv) == rownames(scores))
all(hiv$pRNA.sqrt==scores$pRNA.sqrt)

hiv$S.Score=scores$S.Score
hiv$G2M.Score=scores$G2M.Score
hiv$Phase=scores$Phase
hiv$log.pRNA.sqrt=-log(hiv$pRNA.sqrt)

p=VlnPlot(hiv, features = 'log.pRNA.sqrt', group.by='Phase')
p
ggsave(p, file='output/figure_s10_celltype_anova.pdf', width=4, height=4)
temp=aov(pRNA.sqrt~Phase, hiv@meta.data)
summary(temp)

VlnPlot(hiv, features = 'G2M.Score', group.by='seurat_clusters')
FeaturePlot(hiv, features = 'G2M.Score', reduction = "umap.rna")
cor.test(hiv$pRNA.sqrt, hiv$G2M.Score, method='spearman')

VlnPlot(hiv, features = 'S.Score', group.by='seurat_clusters')
FeaturePlot(hiv, features = 'S.Score', reduction = "umap.rna")
cor.test(hiv$pRNA.sqrt, hiv$S.Score, method='spearman')

plot(hiv$S.Score, hiv$pRNA.sqrt)
smoothScatter(hiv$S.Score, hiv$pRNA.sqrt)

cor.test(hiv$pRNA.sqrt[hiv$pRNA.sqrt>0], hiv$S.Score[hiv$pRNA.sqrt>0], method='spearman')
smoothScatter(hiv$S.Score[hiv$pRNA.sqrt>0], hiv$pRNA.sqrt[hiv$pRNA.sqrt>0])




library(tidyverse)
library(ggplot2)
library(patchwork)
library(Seurat)
library(SeuratObject)
library(Signac)

dir.out <- "derived_data"
dir.fig <- "figures"
file.in <- "hiv_predict.rds"

file.umap <- "umap.pdf"

seurat <- file.path(dir.out, file.in) %>% readRDS()
names(seurat@meta.data)
atac <- seurat@assays$ATAC@counts
tot.atac <- colSums(atac)
d <- tibble(atac = tot.atac, cluster = seurat$seurat_clusters)
mean <- d %>% group_by(cluster) %>% summarise(mean = mean(atac))
mean
# prop.atac <- apply(atac, 2, function(x) sum(x > 0)/nrow(atac))

file.path(dir.fig, file.umap) %>% pdf(width = 10.5, height = 10.5)
p1 <- DimPlot(seurat, reduction = "umap.rna", group.by = "group",
              label = TRUE, label.size = 2.5, repel = TRUE) +
    ggtitle("RNA")
p2 <- DimPlot(seurat, reduction = "umap.atac", group.by = "group",
              label = TRUE, label.size = 2.5, repel = TRUE) +
    ggtitle("ATAC")
p3 <- DimPlot(seurat, reduction = "wnn.umap", group.by = "group",
              label = TRUE, label.size = 2.5, repel = TRUE) +
    ggtitle("WNN")
p4 <- DimPlot(seurat, reduction = "umap.rna", group.by = "seurat_clusters",
              label = TRUE, label.size = 2.5, repel = TRUE) +
    ggtitle("RNA")
p5 <- DimPlot(seurat, reduction = "umap.atac", group.by = "seurat_clusters",
              label = TRUE, label.size = 2.5, repel = TRUE) +
    ggtitle("ATAC")
p6 <- DimPlot(seurat, reduction = "wnn.umap", group.by = "seurat_clusters",
              label = TRUE, label.size = 2.5, repel = TRUE) +
    ggtitle("WNN")
p7 <- FeaturePlot(seurat, features = "pRNA.sqrt", reduction = "umap.rna", order = TRUE)
p8 <- FeaturePlot(seurat, features = "pRNA.sqrt", reduction = "umap.atac", order = TRUE)
p9 <- FeaturePlot(seurat, features = "pRNA.sqrt", reduction = "wnn.umap", order = TRUE)
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + plot_layout(ncol = 3) &
    NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

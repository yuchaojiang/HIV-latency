library(Seurat)
library(SeuratObject)
library(Signac)
# library(EnsDb.Hsapiens.v86)
# library(GenomeInfoDb)
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

dir.out <- "derived_data"
dir.fig <- "figures"
file.in <- "hiv_predict.rds"
file.out <- "no_hiv.rds"
file.umap <- "umap_no_hiv.pdf"
file.vln <- "vln_no_hiv.pdf"

seurat <- file.path(dir.out, file.in) %>% readRDS()

dim(seurat@assays$RNA)
dim(seurat@assays$SCT)
dim(seurat@assays$ATAC)
tail(rownames(seurat@assays$RNA), n = 10)
tail(rownames(seurat@assays$SCT), n = 10)
tail(rownames(seurat@assays$ATAC), n = 10)

feat.rna <- rownames(seurat@assays$RNA)
feat.sct <- rownames(seurat@assays$SCT)
feat.atac <- rownames(seurat@assays$ATAC)
# sapply(list(feat.rna, feat.sct, feat.atac), length)

feat.rna <- feat.rna[-grep("^2D10-", feat.rna)]
# feat.sct <- feat.sct[-grep("^2D10-", feat.sct)]
feat.atac <- feat.atac[-grep("^HIV-", feat.atac)]
# sapply(list(feat.rna, feat.sct, feat.atac), length)

head(feat.rna)

DefaultAssay(seurat) <- "RNA"
seurat
seurat <- seurat[c(feat.rna, feat.atac), ]

names(seurat@meta.data)

# perform RNA analysis
DefaultAssay(seurat) <- "RNA"
seurat <- SCTransform(seurat, verbose = FALSE, return.only.var.genes = FALSE) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
seurat

# perform ATAC analysis
# exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(seurat) <- "ATAC"
seurat <- RunTFIDF(seurat) %>%
    FindTopFeatures(min.cutoff = 'q0') %>%
    RunSVD() %>%
    RunUMAP(reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# perform WNN analysis
seurat <- FindMultiModalNeighbors(seurat, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50)) %>%
    RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
    FindClusters(graph.name = "wsnn", algorithm = 3, verbose = FALSE)

seurat$celltype <- seurat$seurat_clusters
file.path(dir.out, file.out) %>% saveRDS(seurat, file = .)

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

file.path(dir.fig, file.vln) %>% pdf(width = 4, height = 8)
p1 <- VlnPlot(seurat, features = "pRNA.sqrt", group.by = "group")
p2 <- VlnPlot(seurat, features = "pRNA.sqrt", group.by = "seurat_clusters")
p1 + p2 + plot_layout(ncol = 1) & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

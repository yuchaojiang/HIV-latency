setwd("/Users/mwen/Documents/Dissertation/HIV/HIV1/downsample/")

#load the HIV1 dataset
#load the HIV2 dataset
load("../rda/hiv_processed.rda")
load("../../HIV2/hiv_processed.rda")

#Now we want to downsample HIV1 so that the number of cells are the same
set.seed(1)
hiv1 <- hiv[,c(sample(which(hiv$group=="DMSO"), 8317, replace = F), 
                 sample(which(hiv$group=="iBET151"), 7932, replace = F),
               sample(which(hiv$group=="Prostratin"), 8827, replace = F),
               sample(which(hiv$group=="SAHA"), 5755 , replace = F))]

dim(hiv)
table(hiv$group)
hiv=SetIdent(hiv, value = 'orig.ident')

pdf("pdf/one/figure1.pdf", width = 5, height = 5)
VlnPlot(hiv, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

smoothScatter((hiv$nCount_ATAC), (hiv$nCount_RNA), pch=16,
              xlab='Total ATAC counts', ylab='Total RNA count')

dev.off()

# RNA analysis
DefaultAssay(hiv) <- "RNA"
hiv <- SCTransform(hiv, verbose = FALSE, return.only.var.genes = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(hiv) <- "ATAC"
hiv <- RunTFIDF(hiv)
hiv <- FindTopFeatures(hiv, min.cutoff = 'q0')
hiv <- RunSVD(hiv)
hiv <- RunUMAP(hiv, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

hiv <- FindMultiModalNeighbors(hiv, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
hiv <- RunUMAP(hiv, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
hiv <- FindClusters(hiv, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
hiv$celltype=hiv$seurat_clusters



head(Idents(hiv))
hiv=SetIdent(hiv, value = 'group')
head(Idents(hiv))
table(hiv$group)


load("../rda/to_share/transcripts.rda")
transcripts = transcripts[match(rownames(hiv@assays$RNA),transcripts$gene_name)]
all(rownames(hiv@assays$RNA)==transcripts$gene_name)
save(hiv, file='hiv.rda')
save(transcripts, file='transcripts.rda')



# Generate some visualizations
pdf("pdf/one/figure4.pdf", width = 15, height = 15)
p1 <- DimPlot(hiv, reduction = "umap.rna", group.by='seurat_clusters', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") + NoLegend()
p2 <- DimPlot(hiv, reduction = "umap.rna", group.by='group', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") + NoLegend()
p3 <- DimPlot(hiv, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")+ NoLegend()
p4 <- DimPlot(hiv, reduction = "umap.atac", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")+ NoLegend()
p5 <- DimPlot(hiv, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p6 <- DimPlot(hiv, reduction = "wnn.umap", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3 +p4 +p5+p6+plot_layout(ncol=2)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))

dev.off()

library(cowplot)
source("../../HIV2/outs/DotPlot.peak.R")



hiv.genes=transcripts[seqnames(transcripts)=='HIV']$gene_name
pdf("pdf/one/figure3.pdf", width = 15, height = 15)
DotPlot(hiv, features=paste0('sct_', hiv.genes), group.by='group')
#VlnPlot(hiv, features=paste0('sct_', hiv.genes), group.by='group')
dev.off()

hiv.peaks=GRangesToString(hiv@assays$ATAC@ranges[seqnames(hiv@assays$ATAC@ranges)=='HIV'])

pdf("pdf/one/figure5.pdf", width = 15, height = 15)
p1=DotPlot.peak(hiv, features=hiv.peaks[1], group.by='group')
p2=DotPlot.peak(hiv, features=hiv.peaks[2], group.by='group')
#p3=DotPlot.peak(hiv, features=hiv.peaks[3], group.by='group')
#p4=DotPlot.peak(hiv, features=hiv.peaks[4], group.by='group')
p1+p2

dev.off()


# Generate some visualization of the HIV
# https://satijalab.org/signac/articles/visualization.html

pdf("pdf/one/figure6.pdf", width = 15, height = 15)
cov_plot <- CoveragePlot(
  object = hiv,
  region = "HIV-1-11000",
  annotation = FALSE,
  peaks = TRUE
)
cov_plot

cov_plot <- CoveragePlot(
  object = hiv,
  region = "HIV-1-11000",
  annotation = FALSE,
  peaks = TRUE,
  group.by = 'seurat_clusters'
)
cov_plot

cov_plot <- CoveragePlot(
  object = hiv,
  region = "HIV-100-1600",
  annotation = FALSE,
  peaks = TRUE
)
cov_plot

cov_plot <- CoveragePlot(
  object = hiv,
  region = "HIV-1-11000",
  annotation = FALSE,
  peaks = TRUE,
  group.by = 'group'
)
cov_plot

dev.off()



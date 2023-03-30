setwd("/Users/mwen/Documents/Dissertation/HIV/HIV2/")

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(patchwork)


# The 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("outs/filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
rm(inputdata.10x)

dim(rna_counts) #[1] 33539 80000  36602 38758
dim(atac_counts) #163066  80000   150800  38758

# Create Seurat object
hiv <- CreateSeuratObject(counts = rna_counts)
hiv[["percent.mt"]] <- PercentageFeatureSet(hiv, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% c(standardChromosomes(grange.counts),'HIV') # chr1-22, chrX, chrY, HIV
atac_counts <- atac_counts[as.vector(grange.use), ]
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.counts

# Read in gene gtf file
annotations <- rtracklayer::import('../HIV1/hiv1/all.genes.gtf.gz')
seqlevelsStyle(annotations) <- 'UCSC'
annotations=annotations[seqnames(annotations) %in% c(paste0('chr',standardChromosomes(grange.counts)),'HIV'),]
genome(annotations) <- "hg38"
transcripts = annotations[annotations$type =="transcript",]

transcripts.exon=annotations[annotations$type=='exon',]
which(is.na(transcripts.exon$gene_name))
which(seqnames(transcripts.exon)=='HIV') 

transcripts.exon[seqnames(transcripts.exon)=='HIV',]
transcripts.exon[seqnames(transcripts.exon)=='HIV',]$transcript_id
transcripts.exon$gene_name[is.na(transcripts.exon$gene_name)]=transcripts.exon$transcript_id[is.na(transcripts.exon$gene_name)]
transcripts.hiv=transcripts.exon[seqnames(transcripts.exon)=='HIV',]

transcripts.exon$gene_name=gsub('_','-',transcripts.exon$gene_name)

transcripts=transcripts[seqnames(transcripts)!='HIV',]
end(transcripts.hiv[1])=max(end(transcripts.hiv))

transcripts.hiv=transcripts.hiv[1]
transcripts.hiv$gene_id='HIV'
transcripts.hiv$gene_name='HIV'
transcripts.hiv$transcript_id='HIV'
transcripts.hiv$source='HIV'
transcripts=c(transcripts.hiv, transcripts) 


# The scRNA-seq read count only has one HIV gene.
# Merge all HIV genes from the annotation

hiv=hiv[!is.na(match(rownames(hiv@assays$RNA), transcripts$gene_name)),]
any(is.na(match(rownames(hiv@assays$RNA), transcripts$gene_name)))
transcripts=transcripts[match(rownames(hiv@assays$RNA), transcripts$gene_name),]
all(rownames(hiv@assays$RNA)==transcripts$gene_name)

length(transcripts); dim(hiv@assays$RNA)
length(grange.counts); dim(atac_counts)

frag.file <- "outs/atac_fragments.tsv.gz"
frag.hiv <- CreateFragmentObject(
  path = frag.file,
  cells = colnames(atac_counts)
)

# Add additional peaks for HIV that are missed by CellRanger
grange.hiv=GRanges(seqnames='HIV',
                   ranges=IRanges(start=c(251, 521, 871), end=c(500, 850, 1100)))
# Re-calculate the count matrix given the new peaks
atac_counts.hiv <- FeatureMatrix(
  fragments = frag.hiv,
  features = grange.hiv,
  cells = colnames(atac_counts)
)
head(rownames(atac_counts))
tail(rownames(atac_counts))
atac_counts=rbind(atac_counts.hiv, atac_counts)

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  # genome = 'hg38', # this results in error because of the HIV genome; we will add later
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
hiv[["ATAC"]] <- chrom_assay
hiv@assays$ATAC@ranges
rm(chrom_assay)

group=substr(colnames(hiv), nchar(colnames(hiv)), nchar(colnames(hiv)))
hiv$groupid=group
#Donor 1 
hiv$group=factor(c('DMSO','iBET151','Prostratin','SAHA')[as.numeric(group)])
#Donor 2
hiv$group=factor(c('DMSO','SAHA','Prostratin','iBET151')[as.numeric(group)])
table(hiv$group)

smoothScatter((hiv$nCount_ATAC), (hiv$nCount_RNA), pch=16, 
              xlab='Total ATAC counts', ylab='Total RNA count')

VlnPlot(hiv, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

hiv <- subset(
  x = hiv,
  subset = nCount_ATAC < 3e4 &
    nCount_ATAC > 1e3 &
    nCount_RNA < 10000 &
    nCount_RNA > 1000 &
    percent.mt < 20 &
    percent.mt > 5
)
dim(hiv)
table(hiv$group)

VlnPlot(hiv, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

smoothScatter((hiv$nCount_ATAC), (hiv$nCount_RNA), pch=16,
              xlab='Total ATAC counts', ylab='Total RNA count')

hiv.mt=hiv[["percent.mt"]] # Since ChrM is removed. We save this and load it back later.

rna_counts=hiv@assays$RNA@counts
dim(rna_counts); length(transcripts)

atac_counts=hiv@assays$ATAC@counts
dim(atac_counts)
rownames(atac_counts)[1:10]

# At least ~5% of cells expressing the gene or have peak accessibility
# Change to 3000 for donor 1 because it has more cells. 
transcripts=transcripts[rowSums(rna_counts)>1500]
rna_counts=rna_counts[rowSums(rna_counts)>1500,]
atac_counts=atac_counts[rowSums(atac_counts)>1500,]

dim(rna_counts); length(transcripts)
dim(atac_counts)

rna_counts_temp=rna_counts>0
transcripts=transcripts[rowSums(rna_counts_temp)>1500]
rna_counts=rna_counts[rowSums(rna_counts_temp)>1500,]
rm(rna_counts_temp)

atac_counts_temp=atac_counts>0
atac_counts=atac_counts[rowSums(atac_counts_temp)>1500,]
rm(atac_counts_temp)

dim(rna_counts); length(transcripts)
dim(atac_counts)

# Create new object
hiv <- CreateSeuratObject(counts = rna_counts)
hiv[["percent.mt"]] = hiv.mt; rm(hiv.mt)

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  # genome = 'hg38', # this results in error because of the HIV genome; we will add later
  fragments = frag.hiv,
  min.cells = 10,
  annotation = annotations
)
hiv[["ATAC"]] <- chrom_assay
hiv@assays$ATAC@ranges

VlnPlot(hiv, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
dim(hiv@assays$RNA)
dim(hiv@assays$ATAC)
rm(atac_counts); rm(rna_counts); rm(atac_counts.hiv); rm(transcripts.hiv); rm(grange.hiv)
rm(chrom_assay); rm(grange.use)



group=substr(colnames(hiv), nchar(colnames(hiv)), nchar(colnames(hiv)))
hiv$groupid=group
#Donor 1 
hiv$group=factor(c('DMSO','iBET151','Prostratin','Vorinostat')[as.numeric(group)])
#Donor 2
hiv$group=factor(c('DMSO','Vorinostat','Prostratin','iBET151')[as.numeric(group)])
table(hiv$group)


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

group=substr(colnames(hiv), nchar(colnames(hiv)), nchar(colnames(hiv)))
hiv$groupid=group
#hiv$group=factor(c('DMSO','SAHA','Prostratin','iBET151')[as.numeric(group)])
#table(hiv$group)


head(Idents(hiv))
hiv=SetIdent(hiv, value = 'group')
head(Idents(hiv))
table(hiv$group)

save(hiv, file='hiv.rda')
save(transcripts, file='transcripts.rda')



# Generate some visualizations
setwd("~/Dropbox/HIV_primary_cell/")
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(patchwork)
load('hiv.rda')
load('transcripts.rda')
project=substr(colnames(hiv), 1, 1)
hiv$project=project

pdf("pdf/one/figure4.pdf", width = 15, height = 15)
p1 <- DimPlot(hiv, reduction = "umap.rna", group.by='seurat_clusters', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") + NoLegend()
p2 <- DimPlot(hiv, reduction = "umap.rna", group.by='group', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") + NoLegend()
p3 <- DimPlot(hiv, reduction = "umap.rna", group.by='project', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") + NoLegend()

p4 <- DimPlot(hiv, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")+ NoLegend()
p5 <- DimPlot(hiv, reduction = "umap.atac", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")+ NoLegend()
p6 <- DimPlot(hiv, reduction = "umap.atac", group.by = "project", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")+ NoLegend()

p7 <- DimPlot(hiv, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p8 <- DimPlot(hiv, reduction = "wnn.umap", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p9 <- DimPlot(hiv, reduction = "wnn.umap", group.by = "project", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3 +p4 +p5+p6+p7+p8+p9+plot_layout(ncol=3)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))

dev.off()



pdf("pdf/one/figure4.1.new.pdf", width = 20, height = 5)
p1 <- DimPlot(hiv, reduction = "umap.rna", group.by='group', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") 
p2 <- DimPlot(hiv, reduction = "umap.atac", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")
p3 <- DimPlot(hiv, reduction = "umap.rna", group.by='seurat_clusters', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") 
p4 <- DimPlot(hiv, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")

p1 + p2  +p3 +p4 +plot_layout(ncol=4) & theme(plot.title = element_text(hjust = 0.5))

dev.off()


pdf("pdf/one/featureplot.pdf", width = 10, height = 10)
DefaultAssay(hiv) <- "RNA"
FeaturePlot(hiv, features = "pRNA.sqrt")
DefaultAssay(hiv) <- "ATAC"
FeaturePlot(hiv, features = "pATAC.sqrt")
dev.off()



library(cowplot)
source("outs/DotPlot.peak.R")

hiv.genes=transcripts[seqnames(transcripts)=='HIV']$gene_name
pdf("pdf/one/figure3.pdf", width = 5, height = 7)
DotPlot(hiv, features=paste0('sct_', hiv.genes), group.by='group')
#VlnPlot(hiv, features=paste0('sct_', hiv.genes), group.by='group')
dev.off()

hiv.ranges= GRanges(seqnames='HIV', ranges=IRanges(start=1, end=10425))
frag.file <- "outs/atac_fragments.tsv.gz"
DefaultAssay(hiv)='ATAC'
frag.hiv <- CreateFragmentObject(
  path = frag.file,
  cells = colnames(hiv)
)
atac_counts=hiv@assays$ATAC@counts
atac_counts.hiv <- FeatureMatrix(
  fragments = frag.hiv,
  features = hiv.ranges,
  cells = colnames(hiv)
)

atac_counts=rbind(atac_counts.hiv, atac_counts)

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  # genome = 'hg38', # this results in error because of the HIV genome; we will add later
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
hiv[["ATAC"]] <- chrom_assay
hiv@assays$ATAC@ranges
hiv.peaks=GRangesToString(hiv@assays$ATAC@ranges[seqnames(hiv@assays$ATAC@ranges)=='HIV'])


pdf("pdf/one/figure5.pdf", width = 15, height = 8)
p1=DotPlot.peak(hiv, features=hiv.peaks[1], group.by='group')
p2=DotPlot.peak(hiv, features=hiv.peaks[2], group.by='group')
p3=DotPlot.peak(hiv, features=hiv.peaks[3], group.by='group')

#p4=DotPlot.peak(hiv, features=hiv.peaks[4], group.by='group')
p1+p2+p3
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



hiv.genes="HIV"
pdf("pdf/one/figure8.pdf", width = 5, height = 10)
DotPlot(hiv, features=paste0('sct_', hiv.genes), group.by='seurat_clusters')
#VlnPlot(hiv, features=paste0('sct_', hiv.genes), group.by='group')
dev.off()




dim(hiv@assays$RNA)
dim(hiv@assays$SCT)
dim(hiv@assays$ATAC)
tail(rownames(hiv@assays$RNA), n = 10)
tail(rownames(hiv@assays$SCT), n = 10)
tail(rownames(hiv@assays$ATAC), n = 10)

feat.rna <- rownames(hiv@assays$RNA)
feat.sct <- rownames(hiv@assays$SCT)
feat.atac <- rownames(hiv@assays$ATAC)
# sapply(list(feat.rna, feat.sct, feat.atac), length)

feat.rna <- feat.rna[-grep("^HIV$", feat.rna)]
# feat.sct <- feat.sct[-grep("^2D10-", feat.sct)]
feat.atac <- feat.atac[-grep("^HIV-", feat.atac)]
# sapply(list(feat.rna, feat.sct, feat.atac), length)"^apple$"

tail(feat.rna)
head(feat.atac)

DefaultAssay(hiv) <- "RNA"
hiv
hiv <- hiv[c(feat.rna, feat.atac), ]

names(hiv@meta.data)

# perform RNA analysis
DefaultAssay(hiv) <- "RNA"
hiv <- SCTransform(hiv, verbose = FALSE, return.only.var.genes = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
hiv

# perform ATAC analysis
# exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(hiv) <- "ATAC"
hiv <- RunTFIDF(hiv) %>%
  FindTopFeatures(min.cutoff = 'q0') %>%
  RunSVD() %>%
  RunUMAP(reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# perform WNN analysis
hiv <- FindMultiModalNeighbors(hiv, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50)) %>%
  RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
  FindClusters(graph.name = "wsnn", algorithm = 3, verbose = FALSE)

hiv$celltype <- hiv$hiv_clusters
file.path(dir.out, file.out) %>% saveRDS(hiv, file = .)

pdf("figure7.3.pdf",width = 10.5, height = 10.5)
p1 <- DimPlot(hiv, reduction = "umap.rna", group.by = "group",
              label = TRUE, label.size = 2.5, repel = TRUE) +
  ggtitle("RNA")
p2 <- DimPlot(hiv, reduction = "umap.atac", group.by = "group",
              label = TRUE, label.size = 2.5, repel = TRUE) +
  ggtitle("ATAC")
p3 <- DimPlot(hiv, reduction = "wnn.umap", group.by = "group",
              label = TRUE, label.size = 2.5, repel = TRUE) +
  ggtitle("WNN")
p4 <- DimPlot(hiv, reduction = "umap.rna", group.by = "seurat_clusters",
              label = TRUE, label.size = 2.5, repel = TRUE) +
  ggtitle("RNA")
p5 <- DimPlot(hiv, reduction = "umap.atac", group.by = "seurat_clusters",
              label = TRUE, label.size = 2.5, repel = TRUE) +
  ggtitle("ATAC")
p6 <- DimPlot(hiv, reduction = "wnn.umap", group.by = "seurat_clusters",
              label = TRUE, label.size = 2.5, repel = TRUE) +
  ggtitle("WNN")
p7 <- FeaturePlot(hiv, features = "pRNA.sqrt", reduction = "umap.rna", order = TRUE)
p8 <- FeaturePlot(hiv, features = "pRNA.sqrt", reduction = "umap.atac", order = TRUE)
p9 <- FeaturePlot(hiv, features = "pRNA.sqrt", reduction = "wnn.umap", order = TRUE)
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + plot_layout(ncol = 3) &
  NoLegend() & theme(plot.title = element_text(hjust = 0.5))

p1 <- VlnPlot(hiv1, features = "pRNA.sqrt", group.by = "group")
p2 <- VlnPlot(hiv1, features = "pRNA.sqrt", group.by = "seurat_clusters", pt.size = 0)
p1 + p2 + plot_layout(ncol = 1) & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

setwd("~/Dropbox/HIV")

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(patchwork)

# The 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("~/Dropbox/Multiome_Shared/data/HIV/outs/filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
hiv <- CreateSeuratObject(counts = rna_counts)
hiv[["percent.mt"]] <- PercentageFeatureSet(hiv, pattern = "^MT-")


# Read in gene gtf file
annotations <- rtracklayer::import('~/Dropbox/Multiome_Shared/data/HIV/GRCh38_HIV/genes/genes.gtf')
seqlevelsStyle(annotations) <- 'UCSC'
annotations=annotations[seqnames(annotations) %in% c(standardChromosomes(grange.counts),'HIV'),]
genome(annotations) <- "hg38"
transcripts=annotations[annotations$type=='transcript',]
transcripts[seqnames(transcripts)=='HIV',]
transcripts$gene_name[is.na(transcripts$gene_name)]=transcripts$gene_id[is.na(transcripts$gene_name)]
transcripts$gene_name=gsub('_','-',transcripts$gene_name)

hiv=hiv[!is.na(match(rownames(hiv@assays$RNA), transcripts$gene_name)),]
any(is.na(match(rownames(hiv@assays$RNA), transcripts$gene_name)))
transcripts=transcripts[match(rownames(hiv@assays$RNA), transcripts$gene_name),]
all(rownames(hiv@assays$RNA)==transcripts$gene_name)


# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% c(standardChromosomes(grange.counts),'HIV') # chr1-22, chrX, chrY, HIV
atac_counts <- atac_counts[as.vector(grange.use), ]
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.counts
grange.counts[seqnames(grange.counts)=='HIV',]

dim(atac_counts)
length(grange.counts)
table(atac_counts[71335,]) # The second HIV peak is removed below 

# Add additional peaks for HIV that are missed by CellRanger
grange.counts = append(grange.counts, GRanges(seqnames='HIV',
                                             ranges=IRanges(start=c(251, 521, 871), end=c(500, 850, 1100))))
grange.counts = sort(grange.counts)
grange.counts

# Re-calculate the count matrix given the new peaks
frag.file <- "~/Dropbox/Multiome_Shared/data/HIV/outs/atac_fragments.tsv.gz"
frag.hiv <- CreateFragmentObject(
  path = frag.file,
  cells = colnames(atac_counts)
)
atac_counts <- FeatureMatrix(
  fragments = frag.hiv,
  features = grange.counts,
  cells = colnames(atac_counts)
)



chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
hiv[["ATAC"]] <- chrom_assay
hiv@assays$ATAC@ranges

VlnPlot(hiv, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

group=substr(colnames(hiv), nchar(colnames(hiv)), nchar(colnames(hiv)))
hiv$groupid=group
hiv$group=factor(c('DMSO','iBET151','Prostratin','SAHA')[as.numeric(group)])


hiv <- subset(
  x = hiv,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 1e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

VlnPlot(hiv, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()


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

p1 <- DimPlot(hiv, reduction = "umap.rna", group.by='group', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") + NoLegend()
p2 <- DimPlot(hiv, reduction = "umap.atac", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")+ NoLegend()
p3 <- DimPlot(hiv, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p4 <- DimPlot(hiv, reduction = "wnn.umap", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3 +p4 +plot_layout(ncol=2)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))

hiv$celltype=hiv$seurat_clusters


tail(rownames(hiv@assays$RNA@counts))
tail(rownames(hiv@assays$SCT@scale.data))
tail(rownames(hiv@assays$ATAC@counts))

hiv.genes=tail(rownames(hiv@assays$SCT@scale.data))
DotPlot(hiv, features=paste0('sct_', hiv.genes), group.by='group')
#VlnPlot(hiv, features=paste0('sct_', hiv.genes), group.by='group')

hiv.peaks=GRangesToString(hiv@assays$ATAC@ranges[seqnames(hiv@assays$ATAC@ranges)=='HIV'])
#VlnPlot(hiv, features=hiv.peaks, group.by='group', ncol = 2)
DotPlot(hiv, features=hiv.peaks, group.by='group')
library(cowplot)
source('DotPlot.peak.R')
p1=DotPlot.peak(hiv, features=hiv.peaks[1], group.by='group')
p2=DotPlot.peak(hiv, features=hiv.peaks[2], group.by='group')
p3=DotPlot.peak(hiv, features=hiv.peaks[3], group.by='group')
p4=DotPlot.peak(hiv, features=hiv.peaks[4], group.by='group')
p1+p2+p3+p4

Idents(hiv)
hiv=SetIdent(hiv, value = 'group')
Idents(hiv)
table(hiv$group)

save(hiv, file='hiv.rda')
save(transcripts, file='transcripts.rda')



# Generate some additional visualizations
setwd("~/Dropbox/HIV")
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(patchwork)
load('processed_data/hiv.rda')


p1 <- DimPlot(hiv, reduction = "umap.rna", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") + NoLegend()
p2 <- DimPlot(hiv, reduction = "umap.atac",  group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")+ NoLegend()
p3 <- DimPlot(hiv, reduction = "wnn.umap", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv WNN (RNA+ATAC)")

p1+p2+p3

# Generate some visualization of the HIV
# https://satijalab.org/signac/articles/visualization.html

cov_plot <- CoveragePlot(
  object = hiv,
  region = "HIV-1-6000",
  annotation = FALSE,
  peaks = TRUE
)
cov_plot

cov_plot <- CoveragePlot(
  object = hiv,
  region = "HIV-100-1600",
  annotation = FALSE,
  peaks = TRUE
)
p=cov_plot
ggsave(p, file='output/supp_figure4_cov.pdf', width=7, height=3.5)




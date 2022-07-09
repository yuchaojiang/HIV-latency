library(tidyverse)
library(patchwork)
library(Seurat)
library(Signac)
library(BiocParallel)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)

dir.r <- "functions"
tmp <- list.files(dir.r, full.names = TRUE, pattern = ".R$") %>%
    lapply(source)
rm(tmp)

args <- commandArgs(TRUE)
file.h5 <- args[1]
file.gtf <- args[2]
file.frag <- args[3]
file.seurat <- args[4]
file.fig <- args[5]

data  <- file.h5 %>% Read10X_h5()

rna.counts <- data %>% getElement("Gene Expression")
atac.counts <- data %>% getElement("Peaks")

seurat <- CreateSeuratObject(counts = rna.counts)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
annotations <- file.gtf %>% rtracklayer::import()
seqlevelsStyle(annotations) <- "UCSC"

grange.counts <- StringToGRanges(rownames(atac.counts), sep = c(":", "-"))
annotations <- annotations[seqnames(annotations) %in% c(standardChromosomes(grange.counts),"HIV"),]
genome(annotations) <- "hg38"
transcripts <- annotations[annotations$type == "transcript",]
transcripts[seqnames(transcripts) == "HIV",]
transcripts$gene_name[is.na(transcripts$gene_name)] <- transcripts$gene_id[is.na(transcripts$gene_name)]
transcripts$gene_name <- gsub("_", "-", transcripts$gene_name)

seurat <- seurat[!is.na(match(rownames(seurat@assays$RNA), transcripts$gene_name)), ]
any(is.na(match(rownames(seurat@assays$RNA), transcripts$gene_name)))
transcripts <- transcripts[match(rownames(seurat@assays$RNA), transcripts$gene_name),]
all(rownames(seurat@assays$RNA) == transcripts$gene_name)

# add in the ATAC-seq data
# use peaks in standard chromosomes only (chr1-22, chrX, chrY, HIV)
grange.counts <- StringToGRanges(rownames(atac.counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% c(standardChromosomes(grange.counts), "HIV")
atac.counts <- atac.counts[as.vector(grange.use), ]
grange.counts <- StringToGRanges(rownames(atac.counts), sep = c(":", "-"))
grange.counts
grange.counts[seqnames(grange.counts) == "HIV", ]
# remove the second peak
grange.counts <- grange.counts[-length(grange.counts)]

# add additional peaks for HIV that were missed by cellranger
grange.counts <- append(grange.counts, GRanges(seqnames = "HIV",
                                               ranges= IRanges(start = c(251, 521, 851), end = c(500, 850, 1100))))
grange.counts <- sort(grange.counts)
grange.counts

# re-calculate the count matrix given the new peaks
frag <- CreateFragmentObject(
  path = file.frag,
  cells = colnames(atac.counts)
)

atac.counts <- FeatureMatrix(
  fragments = frag,
  features = grange.counts,
  cells = colnames(atac.counts)
)

chrom.assay <- CreateChromatinAssay(
  counts = atac.counts,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = file.frag,
  min.cells = 10,
  annotation = annotations
)

seurat[["ATAC"]] <- chrom.assay

p1 <- VlnPlot(seurat, features = c("nCount_RNA", "nCount_ATAC", "percent.mt"), ncol = 3,
    log = TRUE, pt.size = 0) + NoLegend()

seurat <- subset(
  x = seurat,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 1e3 &
    nCount_RNA < 2.5e4 &
    nCount_RNA > 1e3 &
    percent.mt < 20
)

p2 <- VlnPlot(seurat, features = c("nCount_RNA", "nCount_ATAC", "percent.mt"), ncol = 3,
    log = TRUE, pt.size = 0) + NoLegend()

file.fig %>% pdf(width = 6, height = 9)
p1 / p2
dev.off()

DefaultAssay(seurat) <- "RNA"
seurat <- SCTransform(seurat, verbose = FALSE) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:50, reduction.name = "umap.rna", reduction.key = "rnaUMAP_")

DefaultAssay(seurat) <- "ATAC"
seurat <- RunTFIDF(seurat) %>%
    FindTopFeatures(min.cutoff = 'q0') %>%
    RunSVD() %>%
    RunUMAP(reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# perfrom WNN analysis
seurat <- FindMultiModalNeighbors(seurat, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50)) %>%
    RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
    FindClusters(graph.name = "wsnn", algorithm = 3, verbose = FALSE)

file.seurat %>% saveRDS(seurat, file = .)


library(tidyverse)
library(GenomicRanges)
library(Seurat)
library(SeuratObject)
library(Signac)
# BiocManager::install("plyranges")
# library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(qvalue)

# file.seurat <- "derived_data/hiv_chromvar.rds"
# file.deg <- "source_data/supp_table1_iBET.RNA.markers.down.csv"
# file.tr <- "../derived_data/transcripts.rds"
# file.frag <- "../source_data/atac_fragments.tsv.gz"
# file.motif <- "derived_data/motifs.iBET151.down.csv"

args <- commandArgs(TRUE)
file.seurat <- args[1]
file.deg <- args[2]
file.tr <- args[3]
file.frag <- args[4]
file.motif <- args[5]

seurat <- file.seurat %>% readRDS()
mat.motif <- seurat@assays$ATAC@motifs@data[1:5, 1:5]
head(seurat@assays$ATAC@motifs@motif.names)

# note that the fdr threshold was set to 0.05
deg <- file.deg %>%
    read.csv() %>%
    pull(X) %>%
    grep("^2D10", ., invert = TRUE,  value = TRUE)

# remove unnecessary metadata columns
# remove the HIV genes
gr.tr <- file.tr %>%
    readRDS() %>%
    select(1:8) %>%
    filter(!grepl("^2D10", gene_name))

# get promoter regions
gr.pr <- gr.tr %>% promoters(upstream = 2e3, downstream = 0)

# https://satijalab.org/signac/reference/findmotifs
# https://satijalab.org/signac/reference/featurematrix

DefaultAssay(seurat) <- "ATAC"
# frag <- CreateFragmentObject(
#     path = file.frag,
#     cells = colnames(seurat)
# )
# mat.atac <- FeatureMatrix(
#     fragments = frag,
#     features = gr.pr,
#     cells = colnames(seurat)
# )

# mat.atac[1:5, 1:5]
# dim(mat.atac)
# length(gr.pr)

# d <- gr.pr %>%
#     as.data.frame() %>%
#     mutate(range = GRangesToString(gr.pr)) %>%
#     select(range, gene_name) %>%
#     filter(range %in% rownames(mat.atac))
#
# range.bg <- d %>% pull(range)
# range.deg <- d %>%
#     filter(gene_name %in% deg) %>%
#     pull(range)

# https://satijalab.org/signac/articles/pbmc_vignette.html

# get ATAC peak regions overlapping the promoter regions
gr.peak <- seurat %>%
    GetAssay(assay = "ATAC") %>%
    GetAssayData(slot = "ranges")
gr.deg <- gr.pr %>%
    filter(gene_name %in% deg) %>%
    reduce()
gr.ol <- gr.peak[overlapsAny(gr.peak, gr.deg)]

range.ol <- GRangesToString(gr.ol)

seurat <- RegionStats(
    object = seurat,
    assay = "ATAC",
    genome = BSgenome.Hsapiens.UCSC.hg38
)

motifs <- FindMotifs(
    object = seurat,
    features = range.ol
)
head(motifs, n = 20)

motifs <- motifs %>%
    select(c(1, ncol(motifs), 2:(ncol(motifs) - 1))) %>%
    mutate(qvalue =  signif(qvalue(p = pvalue)$qvalues, 4))
file.motif %>%
    write.csv(motifs, file = ., row.names = FALSE, quote = FALSE)


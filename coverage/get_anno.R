# load packages
# note that dplyr (tidyverse) needs to be loaded first
# https://github.com/lzcyzm/exomePeak/issues/1
library(tidyverse)
# library(Seurat)
# library(Signac)
# library(BiocParallel)
library(GenomicRanges)

dir.src <- "source_data"
dir.drv <- "derived_data"
dir.fig <- "figures"
dir.r <- "functions"

# getGroupedFragment()
# getBedGraph()
# normalizeBedGraph()
# getPlotRegion()

file.in <- "transcripts.rds"
file.out <- "bed.anno.rds"

gr <- file.path(dir.drv, file.in) %>%
    readRDS()
gr <- gr[seqnames(gr) == "HIV", "gene_name"]
bed <- data.frame(gr)
bed$score <- 0
bed$type <- "exon"
bed <- bed[, c("seqnames", "start", "end", "gene_name", "score", "strand", "type")]
names(bed)[4] <- "gene"
bed$gene <- gsub("2D10-", "", bed$gene)

# save the object
file.path(dir.drv, file.out) %>%
    saveRDS(bed, .)

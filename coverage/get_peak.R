# load packages
# note that dplyr (tidyverse) needs to be loaded first
# https://github.com/lzcyzm/exomePeak/issues/1
library(tidyverse)
library(Seurat)
library(Signac)
library(BiocParallel)
library(GenomicRanges)

dir.src <- "source_data"
dir.drv <- "derived_data"
dir.fig <- "figures"
dir.r <- "functions"

# getGroupedFragment()
# getBedGraph()
# normalizeBedGraph()
# getPlotRegion()

file.seurat <- "hiv_predict.rds"
file.bed <- "bed.peaks.rds"

object <- file.path(dir.drv, file.seurat) %>%
    readRDS()

gr.peaks <- object@assays$ATAC@ranges
gr.peaks <- gr.peaks[seqnames(gr.peaks) == "HIV"]
bed.peaks <- data.frame(gr.peaks)
bed.peaks <- bed.peaks[, 1:3]
bed.peaks

file.path(dir.drv, file.bed) %>%
    saveRDS(bed.peaks, .)

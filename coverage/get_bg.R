# note that dplyr (tidyverse) needs to be loaded first
# https://github.com/lzcyzm/exomePeak/issues/1
library(tidyverse)
library(Seurat)
library(Signac)
library(BiocParallel)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(TFBSTools)
library(JASPAR2020)
library(patchwork)
library(nbpMatching)
library(Sushi)
library(grid)
library(gridBase)
library(RColorBrewer)

dir.src <- "source_data"
dir.drv <- "derived_data"
dir.fig <- "figures"
dir.r <- "functions"

# getGroupedFragment()
# getBedGraph()
# normalizeBedGraph()
# getPlotRegion()

# create a gene-specific directory
# dir.g <- file.path(dir.out, gene.name)
# if (!dir.exists(dir.g)) {dir.create(dir.g)}

prefix <- "GRCh38_pHR_"
suffix <- "_gex_HIV_bg.txt"
file.seurat <- "hiv_predict.rds"
file.out <- "list.bg.rds"

list.bg <- list.files(dir.drv, pattern = suffix, full.names = TRUE) %>%
    lapply(read.delim, header = FALSE)
lapply(list.bg, head)
names <- list.files(dir.drv, pattern = suffix, full.names = FALSE) %>%
    gsub(prefix, "", .) %>%
    gsub(suffix, "", .)
names(list.bg) <- names

# split the data frame based on cell types
object <- file.path(dir.drv, file.seurat) %>% readRDS()

# normalize the fragment counts against the number of cells
str(object@meta.data)
class(object$group)
levels(object$group)

d <- object@meta.data %>%
    as_tibble() %>%
    dplyr::select(nCount_RNA, group) %>%
    group_by(group) %>%
    summarise(sum = sum(nCount_RNA))
d

# num.cells <- as.vector(table(object$group))
# names(num.cells) <- names(table(object$celltype))
list.norm <- list()
for (i in 1:length(list.bg)) {
    signal <- list.bg[[i]]
    signal$V4 <- signal$V4/d$sum[i] * 1e6 # signal$V4 <- signal$V4/num.cells[i]
    list.norm[[i]] <- signal
}
rm(signal)
names(list.norm) <- names(list.bg)
lapply(list.norm, head)
# save the object
file.path(dir.drv, file.out) %>%
    file.exists()
file.path(dir.drv, file.out) %>%
    saveRDS(list.norm, .)

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

# create a gene-specific directory
# dir.g <- file.path(dir.out, gene.name)
# if (!dir.exists(dir.g)) {dir.create(dir.g)}

# prefix <- "bg_cluster"
# suffix <- ".txt"
# file.seurat <- "hiv_predict.rds"
# file.seurat <- "seurat.rds"
# file.out <- "list.bg.cluster.rds"
#
# file.exists(file.path(dir.drv, file.out))
# if (file.exists(file.path(dir.drv, file.out))) {
#     stop("The output file already exists.")
# }

args <- commandArgs(TRUE)
file.seurat <- args[1]
file.out <- args[2]
files.bg <- args[3:length(args)]

list.bg <- files.bg %>%
    lapply(read.delim, header = FALSE)
names(list.bg) <- seq_along(list.bg) - 1

object <- file.seurat %>% readRDS()

# normalize the fragment counts against the total read counts
# normalize the fragment counts against the number of cells
str(object@meta.data)
class(object$seurat_clusters)
levels(object$seurat_clusters)

d <- object@meta.data %>%
    as_tibble() %>%
    dplyr::select(nCount_RNA, seurat_clusters) %>%
    group_by(seurat_clusters) %>%
    summarise(sum = sum(nCount_RNA))
d

# num.cells <- as.vector(table(object$group))
# names(num.cells) <- names(table(object$celltype))

# correct the order of the list
ord <- seq_along(list.bg) - 1
ord
list.bg <- list.bg[as.character(ord)]
lapply(list.bg, head)
all(d$seurat_clusters == names(list.bg))

list.norm <- list()
for (i in seq_along(list.bg)) {
    signal <- list.bg[[i]]
    signal$V4 <- signal$V4/d$sum[i] * 1e6 # signal$V4 <- signal$V4/num.cells[i]
    list.norm[[i]] <- signal
}
rm(signal)
names(list.norm) <- names(list.bg)
lapply(list.norm, head)
# save the object
file.out %>%
    saveRDS(list.norm, .)

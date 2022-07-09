library(tidyverse)
library(Seurat)
library(Signac)
library(BiocParallel)
library(GenomicRanges)

# set directories
# dir.src <- "../source_data"
# dir.drv <- "derived_data"
# dir.fig <- "figures"
# dir.r <- "functions"

# set output file names
file.frag <- "derived_data/frag_combined_filtered.txt"
file.seurat <- "derived_data/hiv_predict.rds"
clust <- 0
file.out <- "derived_data/frag_clustter_0.txt"

args <- commandArgs(TRUE)
file.frag <- args[1]
file.seurat <- args[2]
clust <- as.integer(args[3])
file.out <- args[4]

# split the data frame based on clusters
frag <- file.frag %>% read.delim(header = FALSE)
object <- file.seurat %>% readRDS()

map.label <- data.frame(
  barcode = colnames(object@assays$RNA),
  label = as.character(object$celltype)
)
head(map.label)

# remove the number of PCR duplicates
frag <- frag[, -5]
colnames(frag) <- c("chrom", "start", "end", "barcode")
merged <- merge(frag, map.label, by.x = 4, by.y = 1)
head(merged)

# remove barcodes
merged <- merged[, -1]
list.frag <- split(x = merged, f = merged$label)
lapply(list.frag, head)
names(list.frag)

i <- clust + 1
d <- list.frag[[i]]
d <- d[, -4]
cluster <- names(list.frag)[i]
file.out %>%
    write.table(d, file = ., sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


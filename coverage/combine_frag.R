library(tidyverse)
library(Seurat)
library(Signac)
library(BiocParallel)
library(GenomicRanges)

# set directories
dir.src <- "../source_data"
dir.drv <- "derived_data"
# dir.fig <- "figures"
# dir.r <- "functions"

file.seurat <- "hiv_predict.rds"
prefix <- "GRCh38_pHR_"
suffix <- "_gex_HIV_frag.txt"

# set output file names
file.frag <- "frag_combined_filtered.txt"
file.bg <- "list.bg.cluster.rds"

args <- commandArgs(TRUE)
# file.seurat <- args[1]

# split into two scripts

list.files(dir.drv)
# create a directory
# dir.g <- file.path(dir.out, gene.name)
# if (!dir.exists(dir.g)) {dir.create(dir.g)}

list.frag <- list.files(dir.drv, pattern = suffix, full.names = TRUE) %>%
    lapply(read.delim, header = FALSE)
lapply(list.frag, head)
names <- list.files(dir.drv, pattern = suffix, full.names = FALSE) %>%
    gsub(prefix, "", .) %>%
    gsub(suffix, "", .)
names(list.frag) <- names

# look at the mapping quality
lapply(list.frag, function(x) table(x$V5))
# $DMSO
#
#     3   255
#   469 18906
#
# $iBET151
#
#     1     3   255
#     1   601 37758
#
# $Prostrat
#
#     3   255
#   594 46888
#
# $SAHA
#
#     3   255
#   850 75440
#

# change the barcode names
list.frag[[2]] <- list.frag[[2]] %>% mutate(V4 = gsub("-1", "-2", V4))
list.frag[[3]] <- list.frag[[3]] %>% mutate(V4 = gsub("-1", "-3", V4))
list.frag[[4]] <- list.frag[[4]] %>% mutate(V4 = gsub("-1", "-4", V4))

# combine the data frame
frag <- do.call("rbind", list.frag)
filtered <- frag %>% filter(V5 == 255)
dim(frag)
dim(filtered)

file.path(dir.drv, file.frag) %>%
    write.table(filtered, file = ., sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

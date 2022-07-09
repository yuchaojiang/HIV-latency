library(tidyverse)
library(Seurat)
library(Signac)
# library(BiocParallel)

args <- commandArgs(trailingOnly = TRUE)
file.seurat <- args[1] # e.g, "derived_data/seurat.SAHA.rds"
file.out <- args[2] # e.g., "derived_data/hiv.rna.SAHA.rds" # a data frame

# file.seurat <- "derived_data/seurat.SAHA.rds"
# file.out <- "derived_data/hiv.rna.SAHA.rds" # a data frame

seurat <- file.seurat %>% readRDS()

# sum of HIV RNA counts
# sum of RNA counts across genes
# proportions of HIV RNA counts
# sum of HIV SCT counts
# sum of HIV counts across genes
# proportion of HIV SCT counts

mat.rna <- seurat %>% GetAssay("RNA") %>% slot("counts")
hiv.rna <- mat.rna[grep("^2D10-", rownames(mat.rna)), ] %>% apply(2, sum)
tot.rna <- mat.rna %>% apply(2, sum)
prop.rna <- hiv.rna %>% `/`(tot.rna) %>% sqrt()

mat.sct <- seurat %>% GetAssay("SCT") %>% slot("data")
hiv.sct <- mat.sct[grep("^2D10-", rownames(mat.sct)), ] %>% apply(2, sum)
tot.sct <- mat.sct %>% apply(2, sum)
prop.sct <- hiv.sct %>% `/`(tot.sct) %>% sqrt()

# plot(hiv.rna, prop.rna)
# plot(hiv.sct, prop.sct)

d <- data.frame(hiv.rna = hiv.rna,
                tot.rna = tot.rna,
                prop.rna = prop.rna,
                hiv.sct = hiv.sct,
                tot.sct = tot.sct,
                prop.sct = prop.sct)
file.out %>% saveRDS(d, file = .)

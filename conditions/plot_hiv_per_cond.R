library(tidyverse)
library(patchwork)
library(Seurat)
library(Signac)

# dir.src <- "source_data"
# dir.drv <- "derived_data"
# dir.fig <- "figures"
# dir.r <- "functions"
# tmp <- list.files(dir.r, full.names = TRUE, pattern = ".R$") %>%
#     lapply(source)
# rm(tmp)

sample.ids <- c("DMSO", "iBET151", "Prostrat", "SAHA")

file.rna <- "figures/hockeystickplot_hiv_rna_per_cond.pdf"
file.atac <- "figures/hockeystickplot_hiv_atac_per_cond.pdf"
# file.sct <- "figures/hockeystickplot_hiv_sct_per_cond.pdf"
# file.tfidf <- "figures/hockeystickplot_hiv_tfidf_per_cond.pdf"
# file.atac.rna <- "figures/scatterplot_hiv_atac_rna_per_cond.pdf"
# file.tfidf.sct <- "figures/scatterplot_hiv_tfidf_sct_per_cond.pdf"

f <- function(x, dir, assay) paste0("hiv.", assay, ".", x, ".rds") %>%
                     file.path(dir, .) %>% readRDS()
# list.rna <- lapply(X = sample.ids, FUN = f, assay = "rna")
# names(list.rna) <- sample.ids
list.rna <- map(sample.ids, ~f(.x, dir = "derived_data/hiv_rna", assay = "rna"))
names(list.rna) <- sample.ids
# map(list.rna, head)
list.atac <- map(sample.ids, ~f(.x, dir = "derived_data/hiv_atac", assay = "atac"))
names(list.atac) <- sample.ids

f.rank <- function(d, sample.id, val, xlab, ylab){
    plot(1:nrow(d), sort(d[, val], decreasing = FALSE),
         xlab = xlab, ylab = ylab,
         main = sample.id, font.main = 1,
         pch = 20)
}

xlab <- "Rank of HIV RNA proportions"
ylab <- "Transformed HIV RNA proportions"
file.rna %>% pdf()
par(mfrow = c(2, 2))
mapply(f.rank, list.rna, sample.ids, val = "prop.rna", xlab = xlab, ylab = ylab)
dev.off()

xlab <- "Rank of HIV ATAC proportions"
ylab <- "HIV ATAC proportions"
file.atac %>% pdf()
par(mfrow = c(2, 2))
mapply(f.rank, list.atac, sample.ids, val = "prop.atac", xlab = xlab, ylab = ylab)
dev.off()

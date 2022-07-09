library(tidyverse)
library(patchwork)
library(Seurat)
library(Signac)
library(qvalue)
# library(BiocParallel)

# args <- commandArgs(trailingOnly = TRUE)
# file.rna <- args[1] # e.g., "derived_data/hiv.rna.SAHA.rds" # a data frame
# file.seurat <- args[2] # e.g., "derived_data/seurat.SAHA.rds"
# file.out.csv <- args[3] # e.g., "derived_data/linked.gene.SAHA.csv"
# file.fig.hist <- args[4] # e.g., "figures/hist_pval_linked_gene_SAHA.pdf"
# file.fig.scatter <- args[5] # e.g., "figures/scatter_linked_gene_SAHA.pdf"

# list.files("derived_data", pattern = ".csv$")
file.csv <- "derived_data/linked.TF.csv"
list.files("derived_data", pattern = ".rds$")
file.fig.scatter <- "figures/scatter_top_TF.pdf"

list.rna <- vector("list", 5)
list.motif <- vector("list", 5)
list.chromvar <- vector("list", 5)

seurat <- "derived_data/seurat.rds" %>%
    readRDS()
seurat.chromvar <- "derived_data/seurat.chromvar.rds" %>%
    readRDS()
slotNames(seurat.chromvar@assays$ATAC@motifs)
head(seurat.chromvar@assays$ATAC@motifs@motif.names)

list.rna[[1]] <- "derived_data/hiv.rna.rds" %>% readRDS() %>% pull(prop.rna)
list.motif[[1]] <- seurat.chromvar@assays$ATAC@motifs@motif.names
list.chromvar[[1]] <- seurat.chromvar@assays$chromvar@data

list.rna[[2]] <- "derived_data/hiv.rna.DMSO.rds" %>% readRDS() %>% pull(prop.rna)
list.rna[[3]] <- "derived_data/hiv.rna.iBET151.rds" %>% readRDS() %>% pull(prop.rna)
list.rna[[4]] <- "derived_data/hiv.rna.Prostrat.rds" %>% readRDS() %>% pull(prop.rna)
list.rna[[5]] <- "derived_data/hiv.rna.SAHA.rds" %>% readRDS() %>% pull(prop.rna)

list.motif[[2]] <- "derived_data/motif.obj.DMSO.rds" %>%
    readRDS() %>% slot("motif.names")
list.motif[[3]] <- "derived_data/motif.obj.iBET151.rds" %>%
    readRDS() %>% slot("motif.names")
list.motif[[4]] <- "derived_data/motif.obj.Prostrat.rds" %>%
    readRDS() %>% slot("motif.names")
list.motif[[5]] <- "derived_data/motif.obj.SAHA.rds" %>%
    readRDS() %>% slot("motif.names")

list.chromvar[[2]] <- "derived_data/chromvar.DMSO.rds" %>%
    readRDS() %>% slot("data")
list.chromvar[[3]] <- "derived_data/chromvar.iBET151.rds" %>%
    readRDS() %>% slot("data")
list.chromvar[[4]] <- "derived_data/chromvar.Prostrat.rds" %>%
    readRDS() %>% slot("data")
list.chromvar[[5]] <- "derived_data/chromvar.SAHA.rds" %>%
    readRDS() %>% slot("data")

top <- file.csv  %>%
    read.csv() %>%
    pull(TF) %>%
    head(n = 6)
top

sample.id <- c("Aggr.", "DMSO", "iBET151", "Prostrat", "SAHA")

# plot for top linked genes
xlab <- "TF motif score"
ylab <- "Proportion of HIV RNA reads"
# xlab <- "RNA exp."
# ylab <- "sqrt(proportion of HIV RNA reads)"

list.motif[[1]][1:5, 1:5]
list.chromvar[[1]][1:5, 1:5]

file.fig.scatter %>% pdf(width = 10, height = 12)
par(mfrow = c(6, 5), mar = c(3, 3, 3, 1), mgp = c(2, 1, 0))
for (TF.name in top) {
    for (i in 1:length(list.chromvar)) {
        motif.names <- list.motif[[i]]
        motif.name <- names(motif.names)[motif.names == TF.name]
        v.score <- list.chromvar[[i]][motif.name,]
        v.rna <- list.rna[[i]]
        test <- cor.test(v.score, v.rna, method = "spearman")
        est <- test$estimate
        pval <- test$p.value
        # main <- paste0(gene.name, ", ", sample.id[i], "\np-val = ", signif(pval, 4))
        main <- paste0(TF.name, ", ", sample.id[i], "\nest. = ",
                       round(est, 3), ", p-val = ", signif(pval, 3))
        smoothScatter(v.score, v.rna,
                      xlab = xlab, ylab = ylab, cex.lab = 0.9, cex.axis = 0.8,
                      main = main, font.main = 1, cex.main = 0.9)
        abline(lm(v.rna ~ v.score), lty = 2, col = "black")
    }
}
dev.off()

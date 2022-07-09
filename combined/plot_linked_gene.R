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

# file.csv <- "derived_data/linked.gene.csv"
file.fig.scatter <- "figures/scatter_top_gene.pdf"

list.seurat <- vector("list", 5)
list.seurat[[1]] <- "derived_data/seurat.rds" %>%
    readRDS()
list.seurat[[2]] <- "derived_data/seurat.DMSO.rds" %>%
    readRDS()
list.seurat[[3]] <- "derived_data/seurat.iBET151.rds" %>%
    readRDS()
list.seurat[[4]] <- "derived_data/seurat.Prostrat.rds" %>%
    readRDS()
list.seurat[[5]] <- "derived_data/seurat.SAHA.rds" %>%
    readRDS()

# top <- file.csv %>%
#     read.csv() %>%
#     filter(est.spe > 0) %>%
#     pull(gene) %>%
#     grep("^2D10-", x = ., invert = TRUE, value = TRUE) %>%
#     head(n = 6)
top <- c("TSPOAP1", "MALAT1", "NEAT1", "ZFPM2-AS1", "ANXA1", "CELF2")

list.sct <- vector("list", 5)
for (i in 1:length(list.seurat)) {
    list.sct[[i]] <- list.seurat[[i]]@assays$SCT@data
}

list.rna <- vector("list", 5)
list.rna[[1]] <- "derived_data/hiv.rna.rds" %>% readRDS() %>% pull(prop.rna)
list.rna[[2]] <- "derived_data/hiv.rna.DMSO.rds" %>% readRDS() %>% pull(prop.rna)
list.rna[[3]] <- "derived_data/hiv.rna.iBET151.rds" %>% readRDS() %>% pull(prop.rna)
list.rna[[4]] <- "derived_data/hiv.rna.Prostrat.rds" %>% readRDS() %>% pull(prop.rna)
list.rna[[5]] <- "derived_data/hiv.rna.SAHA.rds" %>% readRDS() %>% pull(prop.rna)

sample.id <- c("Aggr.", "DMSO", "iBET151", "Prostrat", "SAHA")

# plot for top linked genes
xlab <- "SCTransform-normalized RNA exp."
ylab <- "Proportion of HIV RNA reads"
# xlab <- "RNA exp."
# ylab <- "sqrt(proportion of HIV RNA reads)"

lapply(list.sct, function(x) "NEAT" %in% rownames(x))

file.fig.scatter %>% pdf(width = 10, height = 12)
par(mfrow = c(6, 5), mar = c(3, 3, 3, 1), mgp = c(2, 1, 0))
for (gene.name in top) {
    print(gene.name)
    for (i in seq_along(list.sct)) {
        v.sct <- list.sct[[i]][gene.name, ]
        v.rna <- list.rna[[i]]
        test <- cor.test(v.sct, v.rna, method = "spearman")
        est <- test$estimate
        pval <- test$p.value
        # main <- paste0(gene.name, ", ", sample.id[i], "\np-val = ", signif(pval, 4))
        main <- paste0(gene.name, ", ", sample.id[i], "\nest. = ",
                       round(est, 3), ", pval = ", signif(pval, 3))
        smoothScatter(v.sct, v.rna,
                      xlab = xlab, ylab = ylab, cex.lab = 0.9, cex.axis = 0.8,
                      main = main, font.main = 1, cex.main = 0.9)
        abline(lm(v.rna ~ v.sct), lty = 2, col = "black")
    }
}
dev.off()

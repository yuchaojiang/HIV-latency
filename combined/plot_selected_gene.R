library(tidyverse)
library(patchwork)
library(Seurat)
library(Signac)
library(qvalue)

file.fig.scatter <- "figures/scatter_selected_gene.pdf"

list.seurat <- vector("list", 5)
list.seurat[[1]] <- "derived_data/hiv_predict.rds" %>%
    readRDS()
list.seurat[[2]] <- "derived_data/seurat.DMSO.rds" %>%
    readRDS()
list.seurat[[3]] <- "derived_data/seurat.iBET151.rds" %>%
    readRDS()
list.seurat[[4]] <- "derived_data/seurat.Prostrat.rds" %>%
    readRDS()
list.seurat[[5]] <- "derived_data/seurat.SAHA.rds" %>%
    readRDS()

top <- c("TSPOAP1", "MALAT1", "NEAT1", "RPL36", "RPS27", "RPL37A")

list.sct <- vector("list", 5)
for (i in 1:length(list.seurat)) {
    list.sct[[i]] <- list.seurat[[i]]@assays$SCT@data
}

list.rna <- vector("list", 5)
list.rna[[1]] <- list.seurat[[1]][[]]$pRNA.sqrt
list.rna[[2]] <- "derived_data/hiv.rna.DMSO.rds" %>% readRDS() %>% pull(prop.rna)
list.rna[[3]] <- "derived_data/hiv.rna.iBET151.rds" %>% readRDS() %>% pull(prop.rna)
list.rna[[4]] <- "derived_data/hiv.rna.Prostrat.rds" %>% readRDS() %>% pull(prop.rna)
list.rna[[5]] <- "derived_data/hiv.rna.SAHA.rds" %>% readRDS() %>% pull(prop.rna)

sample.id <- c("Aggr.", "DMSO", "iBET151", "Prostrat", "SAHA")

# generate plots for selected linked genes
# xlab <- "SCTransform-normalized RNA exp."
# ylab <- "Proportion of HIV RNA reads"
xlab <- "RNA exp."
ylab <- "sqrt(proportion of HIV RNA reads)"

file.fig.scatter %>% pdf(width = 10, height = 12)
par(mfrow = c(6, 5), mar = c(3, 3, 3, 1), mgp = c(2, 1, 0))
j <- 0
for (gene.name in top) {
    j <- j + 1
    for (i in 1:length(list.sct)) {
        # linked.gene <- list.link[[i]]
        v.sct <- list.sct[[i]][gene.name, ]
        v.rna <- list.rna[[i]]
        test <- cor.test(v.sct, v.rna, method = "spearman")
        est <- test$estimate
        # pval <- test$p.value
        fit <- lm(v.rna ~ v.sct)
        pval <- summary(fit)$coefficients[2, 4]
        # main <- paste0(gene.name, ", ", sample.id[i], "\np-val = ", signif(pval, 4))
        main <- paste0(gene.name, ", ", sample.id[i], "\nr = ",
                       round(est, 3), ", p-val = ", signif(pval, 3))
        smoothScatter(v.sct, v.rna,
                      xlab = xlab, ylab = ylab, cex.lab = 0.9, cex.axis = 0.8,
                      main = main, font.main = 1, cex.main = 0.9)
        abline(lm(v.rna ~ v.sct), lty = 2, col = "black")
        at <- min(v.sct) - 0.2 * (max(v.sct) - min(v.sct))
        if (i == 1) {
            mtext(LETTERS[j], at = at, line = 1, font = 1, cex = 1.5)
        }
    }
}
dev.off()

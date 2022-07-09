library(tidyverse)
library(patchwork)
library(Seurat)
library(Signac)
library(qvalue)

file.fig.scatter <- "figures/scatter_selected_TF.pdf"

list.rna <- vector("list", 5)
list.motif <- vector("list", 5)
list.chromvar <- vector("list", 5)

seurat <- "derived_data/hiv_predict.rds" %>%
    readRDS()
seurat.chromvar <- "derived_data/hiv_chromvar.rds" %>%
    readRDS()

list.rna[[1]] <- seurat$pRNA.sqrt
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

top <- c("ELF1", "GABPA", "ETV6", "SNAI2", "SNAI3", "TBX5")
sample.id <- c("Aggr.", "DMSO", "iBET151", "Prostrat", "SAHA")

# plot for top linked genes
xlab <- "TF motif score"
# ylab <- "Proportion of HIV RNA reads"
ylab <- "sqrt(proportion of HIV RNA reads)"

head(list.motif[[1]])
list.chromvar[[1]][1:5, 1:5]

file.fig.scatter %>% pdf(width = 10, height = 12)
par(mfrow = c(6, 5), mar = c(3, 3, 3, 1), mgp = c(2, 1, 0))
j <- 0
for (TF.name in top) {
    j <- j + 1
    for (i in 1:length(list.chromvar)) {
        motif.names <- list.motif[[i]]
        motif.name <- names(motif.names)[motif.names == TF.name]
        v.score <- list.chromvar[[i]][motif.name,]
        v.rna <- list.rna[[i]]
        test <- cor.test(v.score, v.rna, method = "spearman")
        est <- test$estimate
#         pval <- test$p.value
        fit <- lm(v.rna ~ v.score)
        pval <- summary(fit)$coefficients[2, 4]
        main <- paste0(TF.name, ", ", sample.id[i], "\nr = ",
                       round(est, 3), ", p-val = ", signif(pval, 3))
        smoothScatter(v.score, v.rna,
                      xlab = xlab, ylab = ylab, cex.lab = 0.9, cex.axis = 0.8,
                      main = main, font.main = 1, cex.main = 0.9)
        abline(lm(v.rna ~ v.score), lty = 2, col = "black")
        at <- min(v.score) - 0.2 * (max(v.score) - min(v.score))
        if (i == 1) {
            mtext(LETTERS[j], at = at, line = 1, font = 1, cex = 1.5)
        }
    }
}
dev.off()

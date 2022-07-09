library(tidyverse)
library(patchwork)
library(Seurat)
library(Signac)
library(qvalue)
# library(BiocParallel)
# library(GenomicRanges)
# library(EnsDb.Hsapiens.v86)
# library(TFBSTools)
# library(JASPAR2020)
# library(patchwork)

# dir.src <- "source_data"
# dir.drv <- "derived_data"
# dir.fig <- "figures"
# dir.r <- "functions"

# tmp <- list.files(dir.r, full.names = TRUE, pattern = ".R$") %>%
#     lapply(source)
# rm(tmp)

args <- commandArgs(trailingOnly = TRUE)
file.rna <- args[1] # e.g., "derived_data/hiv.rna.SAHA.rds" # a data frame
file.motif <- args[2] # e.g., "derived_data/motif.obj.SAHA.rds"
file.chromvar <- args[3] # e.g., "derived_data/chomvar.SAHA.rds"
file.out.csv <- args[4] # e.g., "derived_data/linked.TF.SAHA.csv"
file.fig.hist <- args[5] # e.g., "figures/hist_pval_linked_TF_SAHA.pdf"
file.fig.scatter <- args[6] # e.g., "figures/scatter_linked_TF_SAHA.pdf"

# file.rna <- "derived_data/hiv.rna.SAHA.rds" # a data frame
# file.motif <- "derived_data/motif.obj.SAHA.rds"
# file.chromvar <- "derived_data/chomvar.SAHA.rds"
# file.out.csv <- "derived_data/linked.TF.SAHA.csv"
# file.fig.hist <- "figures/hist_pval_linked_TF_SAHA.pdf"
# file.fig.scatter <- "figures/scatter_linked_TF_SAHA.pdf"

d.rna <- file.rna %>% readRDS()
mat.motif <- file.chromvar %>% readRDS() %>% slot("data")
all(rownames(d.rna) == colnames(mat.motif))

TF <- file.motif %>%
    readRDS() %>%
    slot("motif.names") %>%
    unlist()
all(names(TF) == rownames(mat.motif))

# note that other quantities can be used
v.rna <- d.rna %>% pull(prop.sct)

# get p-value for each TF
n <- length(TF)

linked.TF <- data.frame(TF = TF,
                        motif = rownames(mat.motif),
                        est.pea = rep(NA, n),
                        pval.pea = rep(NA, n),
                        est.spe = rep(NA, n),
                        pval.spe = rep(NA, n))

for (i in 1:nrow(linked.TF)) {
    test.pea <- cor.test(mat.motif[i, ], v.rna)
    test.spe <- cor.test(mat.motif[i, ], v.rna, method = "spearman")
    linked.TF$est.pea[i] <- test.pea$estimate
    linked.TF$pval.pea[i] <- test.pea$p.value
    linked.TF$est.spe[i] <- test.spe$estimate
    linked.TF$pval.spe[i] <- test.spe$p.value
}

file.fig.hist %>% pdf(width = 8, height = 4)
par(mfrow = c(1, 2))
hist(linked.TF$pval.pea,
     xlab = "Pearson p-value",
     main = "Pearson p-val distribution",
     font.main = 1)
hist(linked.TF$pval.spe,
     xlab = "Spearman p-value",
     main = "Spearman p-val distribution",
     font.main = 1)
dev.off()

linked.TF$qval.pea <- qvalue(p = linked.TF$pval.pea)$qvalues
linked.TF$qval.spe <- qvalue(p = linked.TF$pval.spe)$qvalues
table(linked.TF$qval.pea < 0.1)
table(linked.TF$qval.spe < 0.1)
head(linked.TF[order(linked.TF$qval.spe, decreasing = FALSE), ], n = 20)

# plot(linked.TF$est.spe, -log(linked.TF$pval.spe, 10),
#      pch = 16, cex = 0.6,
#      xlab = "Spearman correlation", ylab = "-log(adj. p-value)",
#      ylim = c(0,20), xlim = c(-0.15, 0.15))
# abline(h = -log(0.05, 10), lty = 2)

linked.TF.ordered <- linked.TF %>% arrange(qval.spe)
linked.TF.ordered %>%
    write.table(file = file.out.csv,
                row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

head(linked.TF.ordered)
top <- linked.TF.ordered %>% slice(1:6) %>% pull(TF)

# plot for top linked TFs
file.fig.scatter %>% pdf(width = 9, height = 6)
par(mfrow = c(2, 3))
for (TF.name in top) {
    pval <- linked.TF$pval.spe[linked.TF$TF == TF.name]
    motif.name <- linked.TF$motif[linked.TF$TF == TF.name]
    smoothScatter(mat.motif[motif.name, ], v.rna,
                  xlab = "TF motif score",
                  ylab = "sqrt(percentage of HIV RNA reads)",
                  main = paste0(TF.name, ", p-val =", signif(pval, 4)), font.main = 1)
    abline(lm(v.rna ~ mat.motif[motif.name, ]), lty = 2, col = "black")
}
par(mfrow = c(1,1))
dev.off()

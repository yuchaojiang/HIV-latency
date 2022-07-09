library(tidyverse)
library(patchwork)
library(Seurat)
library(Signac)
library(qvalue)
# library(BiocParallel)

args <- commandArgs(trailingOnly = TRUE)
file.rna <- args[1] # e.g., "derived_data/hiv.rna.SAHA.rds" # a data frame
file.seurat <- args[2] # e.g., "derived_data/seurat.SAHA.rds"
file.out.csv <- args[3] # e.g., "derived_data/linked.gene.SAHA.csv"
file.fig.hist <- args[4] # e.g., "figures/hist_pval_linked_gene_SAHA.pdf"
file.fig.scatter <- args[5] # e.g., "figures/scatter_linked_gene_SAHA.pdf"

# file.rna <- "derived_data/hiv.rna.SAHA.rds" # a data frame
# file.seurat <- "derived_data/seurat.SAHA.rds"
# file.out.csv <- "derived_data/linked.gene.SAHA.csv"
# file.fig.hist <- "figures/hist_pval_linked_gene_SAHA.pdf"
# file.fig.scatter <- "figures/scatter_linked_gene_SAHA.pdf"

d.rna <- file.rna %>% readRDS()
seurat <- file.seurat %>% readRDS()
DefaultAssay(seurat) <- "SCT"
hvg <- VariableFeatures(seurat)
# remove HIV transcripts
hvg <- hvg %>% grep(pattern = "^2D10-", x = ., value = TRUE, invert = TRUE)

sct <- seurat %>% GetAssay("SCT") %>% slot("data") %>% as.matrix()

# note that other quantities can be used
v.rna <- d.rna %>% pull(prop.sct)

# get p-value for each gene
n <- length(hvg)

linked.gene <- data.frame(gene = hvg,
                          est.pea = rep(NA, n),
                          pval.pea = rep(NA, n),
                          est.spe = rep(NA, n),
                          pval.spe = rep(NA, n))

for (i in 1:nrow(linked.gene)) {
    gene.name <- hvg[i]
    test.pea <- cor.test(sct[gene.name, ], v.rna)
    test.spe <- cor.test(sct[gene.name, ], v.rna, method = "spearman")
    linked.gene$est.pea[i] <- test.pea$estimate
    linked.gene$pval.pea[i] <- test.pea$p.value
    linked.gene$est.spe[i] <- test.spe$estimate
    linked.gene$pval.spe[i] <- test.spe$p.value
}

file.fig.hist %>% pdf(width = 8, height = 4)
par(mfrow = c(1, 2))
hist(linked.gene$pval.pea,
     xlab = "Pearson p-value",
     main = "Pearson p-val distribution",
     font.main = 1)
hist(linked.gene$pval.spe,
     xlab = "Spearman p-value",
     main = "Spearman p-val distribution",
     font.main = 1)
dev.off()
file.fig.hist

linked.gene$qval.pea <- qvalue(p = linked.gene$pval.pea)$qvalues
linked.gene$qval.spe <- qvalue(p = linked.gene$pval.spe)$qvalues
# table(linked.gene$qval.pea < 0.1)
# table(linked.gene$qval.spe < 0.1)
# head(linked.gene[order(linked.gene$qval.spe, decreasing = FALSE), ], n = 10)

linked.gene.ordered <- linked.gene %>% arrange(qval.spe)
linked.gene.ordered %>%
    write.table(file = file.out.csv,
                row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

# head(linked.gene.ordered)
top <- linked.gene.ordered %>% slice(1:6) %>% pull(gene)

# plot for top linked genes
xlab <- "Transformed RNA exp."
ylab <- "Proportion of transformed HIV exp."
# xlab <- "RNA exp."
# ylab <- "sqrt(percentage of HIV RNA reads)"
file.fig.scatter %>% pdf(width = 9, height = 6)
par(mfrow = c(2, 3))
for (gene.name in top) {
    pval <- linked.gene$pval.spe[linked.gene$gene == gene.name]
    smoothScatter(sct[gene.name, ], v.rna,
                  xlab = xlab,
                  ylab = ylab,
                  main = paste0(gene.name, ", p-val =", signif(pval, 4)), font.main = 1)
    abline(lm(v.rna ~ sct[gene.name, ]), lty = 2, col = "black")
}
par(mfrow = c(1,1))
dev.off()

# top <- c("TSPOAP1", "MALAT1", "NEAT1", "ZFPM2-AS1")) {


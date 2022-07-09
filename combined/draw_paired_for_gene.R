library(tidyverse)
library(patchwork)
library(Seurat)
library(Signac)
library(qvalue)

dir.r <- "functions"
tmp <- list.files(dir.r, full.names = TRUE, pattern = ".R$") %>%
    lapply(source)
rm(tmp)

file.fig <- "figures/paired_for_genes.pdf"
file.mat <- "derived_data/mat.transformed.pval.gene.rds"

list.files("derived_data", pattern = "csv", full.name = TRUE) %>%
    grep("linked[.]gene", ., value = TRUE)
list.link <- vector("list", 5)
# list.link[[1]] <- "derived_data/linked.gene.aggr.csv" %>%
#     read.csv()
list.link[[1]] <- "derived_data/linked.gene.csv" %>%
    read.csv()
list.link[[2]] <- "derived_data/linked.gene.DMSO.csv" %>%
    read.csv()
list.link[[3]] <- "derived_data/linked.gene.iBET151.csv" %>%
    read.csv()
list.link[[4]] <- "derived_data/linked.gene.Prostrat.csv" %>%
    read.csv()
list.link[[5]] <- "derived_data/linked.gene.SAHA.csv" %>%
    read.csv()

# take p-values from simple linear regression

lapply(list.link, nrow)
genes <- list.link %>%
    lapply(function(x) x$gene) %>%
    Reduce("intersect", .)
length(genes)

sel <- c("pval", "pval.pea")
tmp <- list.link %>%
    lapply(function(x) x[match(genes, x$gene), ]) %>%
    lapply(function(x) x[, colnames(x) %in% sel]) %>%
    do.call(cbind, .)

cauchy.comb <- function(pvals){
  pvals[pvals > 0.999] <- 0.999
  is.small <- (pvals < 1e-15)
  pvals[!is.small] <- tan((0.5 - pvals[!is.small]) * pi)
  pvals[is.small] <- 1/pvals[is.small]/pi
  cct.stat <- mean(pvals)
  pval <- pcauchy(cct.stat, lower.tail = FALSE)
  return(pval)
}

# combine p-values from the LRA-treatmented samples
pval.comb <- apply(tmp[, -(1:2)], 1, cauchy.comb)
mat.pval <- cbind(tmp, pval.comb)
mat.transformed <- sqrt(-log(mat.pval))
colnames(mat.transformed) <- c("Aggregated", "DMSO", "iBET151", "Prostrat", "SAHA", "Combined")

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor = NULL, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    filter <- !is.na(y) & !is.na(x) & !is.nan(x) & !is.nan(y) & !is.infinite(x) & !is.infinite(y)
    r <- cor(x[filter], y[filter])
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, 'r = ', txt)
    cex.cor <- 1.2 # 0.5/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor) # text(0.5, 0.5, txt, cex = cex.cor * r*2)
}

file.fig %>% pdf()
pairs(
    mat.transformed,
    upper.panel = panel.cor,
    lower.panel = function(x, y){smoothScatter(x, y, add = TRUE)}
)
dev.off()

file.mat %>% saveRDS(mat.transformed, file = .)


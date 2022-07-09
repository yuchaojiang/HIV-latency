# load packages
# note that dplyr (tidyverse) needs to be loaded first
# https://github.com/lzcyzm/exomePeak/issues/1
library(tidyverse)
library(Seurat)
library(Signac)
library(BiocParallel)
library(GenomicRanges)
# library(EnsDb.Hsapiens.v86)
# library(TFBSTools)
# library(JASPAR2020)
# library(patchwork)
# library(nbpMatching)
library(Sushi)
library(grid)
library(gridBase)
# library(RColorBrewer)

dir.src <- "source_data"
dir.drv <- "derived_data"
dir.fig <- "figures"
dir.r <- "functions"

tmp <- list.files(dir.r, full.names = TRUE, pattern = ".R$") %>%
    lapply(source)
rm(tmp)

# getGroupedFragment()
# getBedGraph()
# normalizeBedGraph()
# getPlotRegion()

file.bg <- "list.bg.cluster.rds"
file.anno <- "bed.anno.rds"
file.peak <- "bed.peaks.rds"
file.plot <- "cov_rna_cluster.pdf"

# read in objects
list.files(dir.drv, pattern = ".rds")
list.bg <- file.path(dir.drv, file.bg) %>% readRDS()
bed.anno <- file.path(dir.drv, file.anno) %>% readRDS()
bed.peaks <- file.path(dir.drv, file.peak) %>% readRDS()

# get ymax
lapply(list.bg, head)
ymax <- max(sapply(list.bg, function(x) max(x$V4)))
ymax

# plot
## set the figure size
widths <- 12
num.cov <- length(list.bg)
height.cov <- 0.6 # cm
heights <- c(
    rep(height.cov*1.6, num.cov),
    rep(height.cov*2, 1), # gene annotation
    rep(height.cov*1, 2)) # gene annotation, ATAC peaks, and genome label
cm.to.inch <- 0.393701
tot.width <- sum(widths)*cm.to.inch
tot.height <- sum(heights)*cm.to.inch
c(tot.width, tot.height)

oma <- c(0, 0, 0, 1) # oma <- c(0, 0, 1, 1)

num.col <- 1
mat <- matrix(1:(num.cov + 3), ncol = num.col, byrow = TRUE)
mat
lay <- layout(
    mat,
    widths = widths,
    heights = heights)
layout.show(lay)

# lay.g <- grid.layout(nrow = 2, ncol = 3)
lay.g <- grid.layout(
    nrow = nrow(mat),
    ncol = ncol(mat),
    widths = unit(widths, "null"),
    heights = unit(heights, "null")
)

# set the plot region
chrom <- "HIV"
chromstart <- 1
chromend <- 12541
# chromend <- max(sapply(list.bg, function(x) max(x$V3)))
# [1] 7336
chromend <- max(bed.anno$end) + 500
chromend

# set colors
n <- length(list.bg)
cols <- gg_color_hue(n)
names(list.bg)

# generate a pdf file
# list.files(dir.fig)
# scale <- 0.9
file.path(dir.fig, file.plot) %>%
    pdf(width = tot.width, height = tot.height)
#    pdf(., width = tot.width, height = tot.height * scale) # change the aspect for the main figure
lay <- layout(
    mat,
    widths = widths,
    heights = heights)

par(oma = oma)
par(mar = c(0, 0, 0, 0))

## plot coverage
par(mar = c(0, 6, 0, 0)) # changed
v.offset <- -1.8
h.offset <- (chromend - chromstart)/100
for (i in 1:length(list.bg)) {
    plotBedgraph(list.bg[[i]], chrom, chromstart, chromend, color = cols[i], range = c(0, ymax))
    segments(chromstart, 0, chromend, 0, col = "gray25", lwd = 0.6)
    mtext(paste("Cluster", names(list.bg)[i]), line = v.offset, at = chromstart - h.offset,
          adj = 1, padj = 0.5, cex = 0.6) # changed
}

## plot gene annotations
par(mar = c(0.4, 6, 0.4, 0))
# h.offset.gene <- ((chromend - chromstart)/100)*(nchar(gene.name) + 4.5)
# h.offset.gene <- ((chromend - chromstart)/100)*(nchar(gene.name) + 4)
v.offset.gene <- -1# need to be changed
# label.TSS <- paste0(gene.name, " TSS")
plotGenes(
    bed.anno,
    chrom,
    chromstart,
    chromend,
    col = "navyblue",
    types = bed.anno$type,
    labeltext = FALSE,
    labeloffset = 0.8,
    fonttype = 1,
    fontsize = 0.6,
    # arrowlength = 0.0005,
    maxrows = 50,
    bheight = 0.25,
    lheight = 0.25,
    height = 0.25,
    plotgenetype = "box",
    xpd = TRUE)
v.offsets <- c(1, 3, 2, 1, 2, 1)
margin.gene <- 80
for (i in 1:nrow(bed.anno)) {
    text(bed.anno$end[i] + margin.gene, v.offsets[i], bed.anno$gene[i], adj = 0, cex = 0.8, xpd = TRUE)
}
# if (as.character(strand(TSS.g)) == "+") {
# 	text(start(TSS.g) - h.offset.gene, 1.05, label.TSS, cex = 0.8, xpd = TRUE)
# } else if (as.character(strand(TSS.g)) == "-") {
# 	text(start(TSS.g) + h.offset.gene, 1.05, label.TSS, cex = 0.8, xpd = TRUE)
# }
mtext("Gene", line = v.offset.gene, at = chromstart - h.offset, adj = 1, padj = 0.5, cex = 0.6)

## plot atac peak regions
# par(mar = c(2, 1, 0.5, 1), mgp = c(2, 0.8, 0))
# par(mar = c(0.25, 6, 0.25, 0)) # par(mar = c(0.2, 6, 0.2, 0))
par(mar = c(0.4, 6, 0.4, 0))
v.offset.atac <- -0.6
plotBed(
    bed.peaks,
    chrom = chrom,
    chromstart = chromstart,
    chromend = chromend,
    color = "gray",
    height = 0.2,
    row = "auto" #, row = 1,
    # wiggle = 0.001
)
mtext("ATAC peak", line = v.offset.atac, at = chromstart - h.offset, adj = 1, padj = 0.5, cex = 0.6)

## label genome
par(mgp = c(0.8, 0.1, 0))
v.offset.genome <- -2 # -2.5
tck.genome <- -0.25
labelgenome(
    chrom,
    chromstart,
    chromend,
    n = 5,
    scale = "bp",
    scalefont = 1,
    scaleline = 10, # scaleline = 0.01,
    scalecex = 0,
    chromfont = 1,
    chromcex = 0,
    chromline = 10, # chromline = 0.01,
    line = 0.4,
    cex.axis = 0.8,
    lwd = 0.6,
    tck = tck.genome)
mtext(chrom, line = v.offset.genome, at = chromstart - h.offset, adj = 1, padj = 0.5, cex = 0.5)
mtext("bp", line = v.offset.genome, at = chromend + h.offset, adj = 1, padj = 0.5, cex = 0.5)

dev.off()

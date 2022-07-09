library(tidyverse)
library(patchwork)
# library(BiocParallel)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

# dir.src <- "source_data"
# dir.drv <- "derived_data"
# dir.fig <- "figures"
# dir.r <- "functions"

# tmp <- list.files(dir.r, full.names = TRUE, pattern = ".R$") %>%
#     lapply(source)
# rm(tmp)

# use the following code in snakefile
# rule run_chromvar:
#     input:
#         "derived_data/seurat.{sample}.rds"
#     output:
#         motif = "derived_data/motif.obj.{sample}.rds",
#         chromvar = "derived_data/chromvar.{sample}.rds"
#     shell:
#         "Rscript run_chromvar.R {input} {output.motif} {output.chromvar}"

args <- commandArgs(trailingOnly = TRUE)
file.seurat <- args[1] # e.g, "derived_data/seurat.SAHA.rds"
file.out.motif <- args[2] # e.g., "derived_data/motif.obj.SAHA.rds"
file.out.chromvar <- args[3] # e.g., "derived_data/chomvar.SAHA.rds"

# file.seurat <- "derived_data/seurat.SAHA.rds"
# file.out <- "derived_data/motif.obj.SAHA.rds"
# file.out.chromvar <- "derived_data/chomvar.SAHA.rds"

seurat <- file.seurat %>% readRDS()

# remove chr HIV so that we can apply chromVAR by its default
seurat@assays$ATAC <- subset(seurat@assays$ATAC,
                             rownames(seurat)[as.character(seqnames(granges(seurat))) != "HIV"])

# scan the DNA sequence of each peak for the presence of each motif, and create a motif object
DefaultAssay(seurat) <- "ATAC"

# note that this process does not depend on the multiome data
# and that it should eventually be in a separate script to avoid repetition
pwm.set <- getMatrixSet(x = JASPAR2020,
                        opts = list(species = 9606, all_versions = FALSE))
gr.motif <- matchMotifs(
    pwms = pwm.set,
    subject = granges(seurat),
    out = "positions",
    genome = "hg38"
)

# note that this is specific to the data
# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
# https://satijalab.org/signac/reference/createmotifmatrix
mat.motif <- CreateMotifMatrix(features = granges(seurat),
                               pwm = pwm.set,
                               genome = "hg38",
                               use.counts = FALSE)
motif.object <- CreateMotifObject(data = mat.motif,
                                  pwm = pwm.set,
                                  positions = gr.motif)
seurat <- SetAssayData(seurat,
                       assay = "ATAC",
                       slot = "motifs",
                       new.data = motif.object)

motif.object %>% saveRDS(file = file.out.motif)

# run chromvar
# note that this step adds a chomvar assay object to the seurat object
# without affecting the motif object
seurat <- RunChromVAR(
    object = seurat,
    genome = BSgenome.Hsapiens.UCSC.hg38
)
seurat %>% GetAssay("chromvar") %>% saveRDS(file = file.out.chromvar)


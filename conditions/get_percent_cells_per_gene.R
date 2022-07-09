library(tidyverse)
library(Seurat)
library(Signac)

# dir.src <- "source_data"
dir.drv <- "derived_data/seurat"
# dir.fig <- "figures"
# dir.r <- "functions"
# tmp <- list.files(dir.r, full.names = TRUE, pattern = ".R$") %>%
#     lapply(source)
# rm(tmp)

sample.ids <- c("DMSO", "iBET151", "Prostrat", "SAHA")

# args <- commandArgs(TRUE)
# file.csv <- args[1]
# files <- args[2:length(args)]

file.csv <- "derived_data/percent.cells.per.gene.csv"
file.range <- "derived_data/percent.cells.per.gene.range.csv"

f <- function(x, dir) paste0("seurat.", x, ".rds") %>%
                     file.path(dir, .) %>% readRDS()

list.seurat <- map(sample.ids, ~f(.x, dir = dir.drv))
names(list.seurat) <- sample.ids

# seurat <- list.seurat[[1]]
# mat.rna <- seurat %>% GetAssay("RNA") %>% slot("counts")
# mat.rna.hiv <- mat.rna[grep("^2D10-", rownames(mat.rna)), ]
# mat.rna.hiv[, 1:5]

f.pc <- function(seurat) {
    mat.rna <- seurat %>% GetAssay("RNA") %>% slot("counts")
    mat.rna[grep("^2D10-", rownames(mat.rna)), ] %>%
        apply(1, function(x) sum(x > 0)) %>%
        `/`(ncol(mat.rna)) %>%
        `*`(100)
}

pc <- map(list.seurat, f.pc) %>%
    do.call("rbind", .)
rownames(pc) <- sample.ids

file.csv
file.csv %>%
    write.csv(pc, file = ., quote = FALSE)

range.pc <- rbind(
    range(pc[1, ]),
    range(pc[2:4, ])
)
colnames(range.pc) <- c("min", "max")
rownames(range.pc) <- c("control", "treatment")
file.range %>%
    write.csv(range.pc, file = ., quote = FALSE)


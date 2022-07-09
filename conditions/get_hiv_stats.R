library(tidyverse)
library(Seurat)
library(Signac)

# dir.src <- "source_data"
# dir.drv <- "derived_data"
# dir.fig <- "figures"
# dir.r <- "functions"
# tmp <- list.files(dir.r, full.names = TRUE, pattern = ".R$") %>%
#     lapply(source)
# rm(tmp)

sample.ids <- c("DMSO", "iBET151", "Prostrat", "SAHA")

# args <- commandArgs(TRUE)
# file.csv <- args[1]
# files <- args[2:length(args)]

file.csv <- "derived_data/percent.exp.range.csv"
# files <- list.files("derived_data", pattern = "^hiv.rna.", full.names = TRUE)
# list.rna <- map(files, readRDS)
# names(list.rna) <- sample.ids

f <- function(x, dir, assay) paste0("hiv.", assay, ".", x, ".rds") %>%
                     file.path(dir, .) %>% readRDS()

list.rna <- map(sample.ids, ~f(.x, dir = "derived_data/hiv_rna/", assay = "rna"))
names(list.rna) <- sample.ids
map(list.rna, head)

f.pc <- function(d) {
    d %>%
        filter(hiv.rna > 0) %>%
        pull(prop.rna) %>%
        range() %>%
        `*`(100)
}

pc <- map(list.rna, f.pc) %>%
    as.data.frame()

file.csv %>%
    write.csv(pc, file = ., quote = FALSE, row.names = FALSE)

library(tidyverse)
library(Seurat)

# dir.src <- "../source_data"
# dir.drv <- "derived_data"
# dir.fig <- "figures"

args <- commandArgs(TRUE)
file.seurat  <- args[1]
group <- args[2]
file.h5 <- args[3]
file.bc <- args[4]

# groups <- c("DMSO", "iBET151", "Prostrat", "SAHA")
# file.seurat <- "derived_data/hiv_predict.rds"
# group <- "DMSO"
# file.h5 <- "../source_data/GRCh38_pHR_DMSO_filtered_feature_bc_matrix.h5"
# file.bc <- paste0("derived_data/bc_group_", group, ".csv")

# get barcodes from the aggregated data
seurat <- file.seurat %>% readRDS()
bc.aggr <- colnames(seurat@assays$RNA)
str(seurat[[]])

# get barcodes from a condition-specific data
bc <- file.h5 %>%
    Read10X_h5 %>%
    getElement(name = "Gene Expression") %>%
    colnames()
head(bc)

if (group == "DMSO") {
    bc <- bc
} else if (group == "iBET151") {
    bc <- bc %>% gsub("-1", "-2", .)
} else if (group == "Prostrat") {
    bc <- bc %>% gsub("-1", "-3", .)
} else if (group == "SAHA") {
    bc <- bc %>% gsub("-1", "-4", .)
} else {
    stop("The group argument is incorrect.")
}

# select barcodes that are also present in the aggregated data
table(!is.na(match(bc, bc.aggr)))
bc <- bc[bc %in% bc.aggr]

if (group == "DMSO") {
    bc <- bc
} else if (group == "iBET151") {
    bc <- bc %>% gsub("-2", "-1", .)
} else if (group == "Prostrat") {
    bc <- bc %>% gsub("-3", "-1", .)
} else if (group == "SAHA") {
    bc <- bc %>% gsub("-4", "-1", .)
}

file.bc %>%
    write.table(bc, file = .,
                row.names = FALSE, col.names = FALSE,
                quote = FALSE, sep = ",")



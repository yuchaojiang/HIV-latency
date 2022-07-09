library(tidyverse)
library(patchwork)
library(Seurat)
library(SeuratObject)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
# package.version("GenomeInfoDb")

args <- commandArgs(trailingOnly = TRUE)
file.h5 <- args[1] # e.g., "source_data/GRCh38_pHR_SAHA_filtered_feature_bc_matrix.h5"
file.gtf <- args[2] # e.g., "source_data/genes.gtf"
file.frag <- args[3] # e.g., "source_data/GRCh38_pHR_SAHA_atac_fragments.tsv.gz"
file.seurat <- args[4] # e.g, "derived_data/seurat.SAHA.rds"
file.out.df <- args[5] # e.g., "derived_data/hiv.atac.SAHA.rds" # a data frame
file.out.seurat <- args[6] # e.g., "derived_data/seurat.hiv.atac.SAHA.rds"

# file.h5 <- "source_data/GRCh38_pHR_SAHA_filtered_feature_bc_matrix.h5"
# file.gtf <- "source_data/genes.gtf"
# file.frag <- "source_data/GRCh38_pHR_SAHA_atac_fragments.tsv.gz"
# file.seurat <- "derived_data/seurat.SAHA.rds"
# file.out.df <- "derived_data/hiv.atac.SAHA.rds" # a data frame
# file.out.seurat <- "derived_data/seurat.hiv.atac.SAHA.rds"

inputdata.10x <- file.h5 %>% Read10X_h5()

# extract RNA and ATAC data
rna.counts <- inputdata.10x$`Gene Expression`
atac.counts <- inputdata.10x$Peaks

# create a Seurat object
seurat <- CreateSeuratObject(counts = rna.counts)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# add in the ATAC-seq data
# use peaks in chr1-22, chrX, chrY, SEURAT
gr.counts <- StringToGRanges(rownames(atac.counts), sep = c(":", "-"))
gr.use <- seqnames(gr.counts) %in% c(standardChromosomes(gr.counts), "HIV")
atac.counts <- atac.counts[as.vector(gr.use), ]
gr.counts <- StringToGRanges(rownames(atac.counts), sep = c(":", "-"))
gr.counts
# grange.counts[seqnames(grange.counts) == "HIV", ]
# dim(atac.counts)
# length(grange.counts)
# check the first HIV peak
# table(atac.counts[length(grange.counts) - 1,])
# check the second HIV peak
# table(atac.counts[length(grange.counts),])

# add additional peaks for HIV that are missed by CellRanger
# chrNameLength.txt
# HIV	12541
# start <- c(251, 521, 871, 1)
# end <- c(500, 850, 1100, 12541)
# add the entire HIV genome region instead
# note that ranges should be non-overlapping for CreateFragmentObject()
# remove the peaks
gr.counts <- gr.counts[-which(seqnames(gr.counts) == "HIV")]
start <- 1
end <- 12541
gr.counts <- append(gr.counts, GRanges(seqnames = "HIV",
                                       ranges = IRanges(start = start, end = end)))
gr.counts <- sort(gr.counts)
gr.counts

# read in gene gtf file
annotations <- file.gtf %>% rtracklayer::import()
seqlevelsStyle(annotations) <- "UCSC"

annotations <- annotations[seqnames(annotations) %in% c(standardChromosomes(gr.counts), ""), ]
genome(annotations) <- "hg38"
transcripts <- annotations[annotations$type == "transcript", ]
transcripts[seqnames(transcripts) == "", ]
transcripts$gene_name[is.na(transcripts$gene_name)] <- transcripts$gene_id[is.na(transcripts$gene_name)]
transcripts$gene_name <- gsub("_", "-", transcripts$gene_name)

seurat <- seurat[!is.na(match(rownames(seurat@assays$RNA), transcripts$gene_name)), ]
any(is.na(match(rownames(seurat@assays$RNA), transcripts$gene_name)))
transcripts <- transcripts[match(rownames(seurat@assays$RNA), transcripts$gene_name), ]
all(rownames(seurat@assays$RNA) == transcripts$gene_name)

# re-calculate the count matrix given the new peaks
frag.seurat <- CreateFragmentObject(
  path = file.frag,
  cells = colnames(atac.counts)
)
atac.counts <- FeatureMatrix(
  fragments = frag.seurat,
  features = gr.counts,
  cells = colnames(atac.counts)
)

# note that the following code requires GenomInfoDb 1.30.0
chrom.assay <- CreateChromatinAssay(
  counts = atac.counts,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = file.frag,
  min.cells = 10,
  annotation = annotations
)

seurat[["ATAC"]] <- chrom.assay

# DefaultAssay(seurat) <- "RNA"
# seurat <- SCTransform(seurat, verbose = FALSE) %>%
#     RunPCA() %>%
#     RunUMAP(dims = 1:50, reduction.name = "umap.rna", reduction.key = "rnaUMAP_")

DefaultAssay(seurat) <- "ATAC"
seurat <- RunTFIDF(seurat) %>%
    FindTopFeatures(min.cutoff = 'q0') %>%
    RunSVD() %>%
    RunUMAP(reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# perfrom WNN analysis
# seurat <- FindMultiModalNeighbors(seurat, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50)) %>%
#     RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
#     FindClusters(graph.name = "wsnn", algorithm = 3, verbose = FALSE)

seurat2 <- file.seurat %>% readRDS()
# all(!is.na(match(colnames(seurat2), colnames(atac.counts))))
seurat <- seurat[, match(colnames(seurat2), colnames(seurat))]
file.out.seurat %>% saveRDS(seurat, file = .)

# sum of HIV ATAC counts
# sum of ATAC counts across genes
# proportions of ATAC RNA counts
# sum of HIV TF-IDF scores
# sum of HIV TF-IDF scores across genes
# proportion of HIV TF-IDF scores

mat.atac <- seurat %>% GetAssay("ATAC") %>% slot("counts")
hiv.atac <- mat.atac[grep("^HIV-", rownames(mat.atac)), , drop = FALSE] %>% apply(2, sum)
tot.atac <- mat.atac %>% apply(2, sum)
prop.atac <- hiv.atac %>% `/`(tot.atac) %>% sqrt()

mat.tfidf <- seurat %>% GetAssay("ATAC") %>% slot("data")
hiv.tfidf <- mat.tfidf[grep("^HIV-", rownames(mat.tfidf)), , drop = FALSE] %>% apply(2, sum)
tot.tfidf <- mat.tfidf %>% apply(2, sum)
prop.tfidf <- hiv.tfidf %>% `/`(tot.tfidf) %>% sqrt()

# plot(hiv.atac, prop.atac)
# plot(hiv.tfidf, prop.tfidf)

d <- data.frame(hiv.atac = hiv.atac,
                tot.atac = tot.atac,
                prop.atac = prop.atac,
                hiv.tfidf = hiv.tfidf,
                tot.tfidf = tot.tfidf,
                prop.tfidf = prop.tfidf)
file.out.df %>% saveRDS(d, file = .)

# renv::snapshot()

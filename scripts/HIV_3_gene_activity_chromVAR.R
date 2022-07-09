setwd("~/Dropbox/HIV")

###############################################
### Get ChromVAR 
###############################################

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

library(glmnet)
library(pheatmap)
library(fields)
library(qvalue)
library(gplots)
library(patchwork)
library(olsrr)
library(pdftools)
library(ape)
library(dendextend)


load('hiv.rda')
dim(hiv@assays$ATAC)
load('transcripts.rda')

# Remove chr HIV so that we can apply chromVAR by its default
hiv@assays$ATAC=subset(hiv@assays$ATAC, rownames(hiv)[as.character(seqnames(granges(hiv)))!='HIV'])

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(hiv) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.positions <- matchMotifs(
  pwms = pwm_set,
  subject = granges(hiv),
  out = 'positions',
  genome = 'hg38'
)

motif.matrix <- CreateMotifMatrix(features = granges(hiv), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set, positions = motif.positions)
motif.object # A Motif object containing 633 motifs in 105913 regions

hiv <- SetAssayData(hiv, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes 
hiv <- RunChromVAR(
  object = hiv,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

dim(hiv@assays$chromvar)

get_motif_umap=function(TF,...){
  motif.name <- ConvertMotifID(hiv, name = TF)
  motif_plot <- FeaturePlot(hiv, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"))
  motif_plot <- motif_plot + ggtitle(paste('chromVAR access. deviation: TF',TF,'motif', names(which(hiv@assays$ATAC@motifs@motif.names==TF))))
  motif_plot
}
get_motif_umap('JUN') | get_motif_umap('FOS')

hiv@assays$chromvar@data[1:5, 1:5] # chromvar data

# Visualization by UMAP
library(reticulate)
use_python('/Users/yuchaojiang/pythonProject/bin/python', required = TRUE)

k=2; distance = "euclidean"; n_neighbors = 10; min_dist = 0.1; rand.seed = 42
set.seed(rand.seed)
reticulate::py_set_seed(rand.seed)
UMAP <- reticulate::import("umap")
umapper <- UMAP$UMAP(n_components = as.integer(k), metric = distance, 
                     n_neighbors = as.integer(n_neighbors), min_dist = min_dist)
Rumap <- umapper$fit_transform
umap <- Rumap(t(as.matrix(hiv@assays$chromvar@data)))
componentX = umap[, 1]
componentY = umap[, 2]

toplot=data.frame(chromVARUMAP1=componentX, chromVARUMAP2=componentY,
                  group=hiv$group)
sp<-ggplot(toplot, aes(x=chromVARUMAP1, y=chromVARUMAP2, color=group)) + geom_point(size = 0.5)+
  guides(colour = guide_legend(override.aes = list(size=5)))
sp
ggsave(sp, file='output/supp_figure1_chromVAR.pdf', width=5, height=4)


VlnPlot(hiv, features=ConvertMotifID(hiv, name = 'JUN')) | VlnPlot(hiv, features=ConvertMotifID(hiv, name = 'FOS'))

save(hiv, file='processed_data/hiv_chromvar.rda')




###############################################
### Get gene activity
###############################################
setwd("~/Dropbox/HIV")

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

library(glmnet)
library(pheatmap)
library(fields)
library(qvalue)
library(gplots)
library(patchwork)
library(olsrr)
library(pdftools)
library(ape)
library(dendextend)


load('processed_data/hiv.rda')
dim(hiv@assays$ATAC)
load('processed_data/transcripts.rda')
transcripts[seqnames(transcripts)=='HIV',]


temp=match(rownames(hiv@assays$RNA),transcripts$gene_name)
temp=temp[!is.na(temp)]
transcripts=transcripts[temp,]
genebodyandpromoter.coords <- Extend(x = transcripts, upstream = 2000, downstream = 200)
rm(temp)

# create a gene by cell matrix
frag.file <- "~/Dropbox/Multiome_Shared/data/HIV/outs/atac_fragments.tsv.gz"
frag.hiv <- CreateFragmentObject(
  path = frag.file,
  cells = colnames(hiv)
)
gene.activities <- FeatureMatrix(
  fragments = frag.hiv,
  features = genebodyandpromoter.coords,
  cells = colnames(hiv)
)
# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- paste0('GeneAct.',gene.key[rownames(gene.activities)])

hiv[["GeneActivity"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(hiv) <- "GeneActivity"

hiv <- SCTransform(hiv,assay="GeneActivity", verbose = FALSE) %>% 
  RunPCA(verbose=F) %>% FindNeighbors(reduction = 'pca', dims = 1:50)%>%
  FindClusters(verbose = FALSE, algorithm = 3)%>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.geneactivity', reduction.key = 'geneactivityUMAP_') 

DimPlot(hiv, reduction = "umap.geneactivity", group.by='group')

# # This is not precise for the HIV genes since 2K upstream is too large for the HIV chromosome
# DotPlot(hiv, features=paste0('sct_GeneAct.',transcripts[seqnames(transcripts)=='HIV']$gene_name), group.by = 'group', scale='FALSE')

toplot=data.frame(gene.activityUMAP1=hiv@reductions$umap.geneactivity@cell.embeddings[,1],
                  gene.activityUMAP2=hiv@reductions$umap.geneactivity@cell.embeddings[,2],
                  group=hiv$group)
sp<-ggplot(toplot, aes(x=gene.activityUMAP1, y=gene.activityUMAP2, color=group)) + geom_point(size = 0.5)+
  guides(colour = guide_legend(override.aes = list(size=5)))
sp

ggsave(sp, file='output/supp_figure1_gene.activity.pdf', width=5, height=4)



#!/usr/bin/env Rscript



###############################################
### Get ChromVAR 
###############################################

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(GenomicRanges)

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


library(BiocParallel) # parallelization
register(MulticoreParam(8)) 


#load('/Users/mwen/Documents/Dissertation/HIV/full/hiv.rda')
#dim(hiv@assays$ATAC)
#[1] 31470 63674

#load('/Users/mwen/Documents/Dissertation/HIV/full/transcripts.rda')

load('../full/hiv.rda')
dim(hiv@assays$ATAC)
load('../full/transcripts.rda')

# Remove chr HIV so that we can apply chromVAR by its default
hiv@assays$ATAC=subset(hiv@assays$ATAC, rownames(hiv)[as.character(seqnames(granges(hiv)))!='HIV'])
hiv@assays$ATAC=subset(hiv@assays$ATAC, rownames(hiv)[as.character(seqnames(granges(hiv)))!='MT'])

temp=granges(hiv1)

#length = 30005 removed 8 
seqlevelsStyle(temp) <- 'UCSC' # Need to add "chr"

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(hiv) <- "ATAC"
# fetches matrix data for all matrices in the database matching criteria defined by the named arguments 
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
#find motif matches 
motif.positions <- matchMotifs(
  pwms = pwm_set,
  subject = temp,
  out = 'positions',
  genome = 'hg38'
)

motif.matrix <- CreateMotifMatrix(features = temp, pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.matrix = motif.matrix[!duplicated(rownames(motif.matrix)),]
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set, positions = motif.positions)
motif.object # A Motif object containing 633 motifs in  29288 regions



# Now need to remove "chr" otherwise the program complains
#rownames(motif.object@data)=gsub('chr','',rownames(motif.object@data))
#seqlevelsStyle(motif.object@positions)='NCBI'
hiv <- SetAssayData(hiv, assay = 'ATAC', slot = 'motifs', new.data = motif.object)
#seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38)='NCBI'


hiv <- RunChromVAR(
  object = hiv,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

dim(hiv@assays$chromvar)

save(hiv, file='hiv_chromvar.rda')




get_motif_umap=function(TF,...){
  motif.name <- ConvertMotifID(hiv, name = TF)
  motif_plot <- FeaturePlot(hiv, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"))
  motif_plot <- motif_plot + ggtitle(paste('chromVAR access. deviation:\n TF',TF,'motif', names(which(hiv@assays$ATAC@motifs@motif.names==TF))))
  motif_plot
}

pdf("pdf/two/figure1.pdf", width = 20, height = 5)
get_motif_umap('JUN(var.2)') | get_motif_umap('FOS') | DimPlot(hiv, reduction = "umap.atac",  group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC") |DimPlot(hiv, reduction = "umap.atac",  group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")
get_motif_umap('JUN(var.2)') | get_motif_umap('FOS') | DimPlot(hiv, reduction = "umap.rna",  group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC") |DimPlot(hiv, reduction = "umap.rna",  group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")
dev.off()

# Visualization by UMAP
library(reticulate)
use_python('/opt/anaconda3/bin/python3.9', required = TRUE)

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

pdf("pdf/two/figure2.pdf", width = 10, height = 10)
VlnPlot(hiv, features=ConvertMotifID(hiv, name = 'JUN(var.2)')) | VlnPlot(hiv, features=ConvertMotifID(hiv, name = 'FOS'))
sp
dev.off()

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
transcripts[seqnames(transcripts)=='HIV',]


temp=match(rownames(hiv@assays$RNA),transcripts$gene_name)
temp=temp[!is.na(temp)]
transcripts=transcripts[temp,]
genebodyandpromoter.coords <- Extend(x = transcripts, upstream = 2000, downstream = 200)
rm(temp)


# create a gene by cell matrix
frag.file <- "outs/atac_fragments.tsv.gz"
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

DimPlot(hiv, reduction = "umap.geneactivity", group.by='group')+ggtitle('')

# This is not precise for the HIV genes since 2K upstream is too large for the HIV chromosome
DotPlot(hiv, features=paste0('sct_GeneAct.',transcripts[seqnames(transcripts)=='HIV']$gene_name), group.by = 'group', scale='FALSE')

save(hiv, file='hiv.gene.activity.rda')


######PLOTS----

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
library(DescTools)
library(destin)



# Visualization by UMAP
library(reticulate)
load("hiv_chromvar.rda")

atac=hiv@assays$ATAC@counts
dim(atac)
rownames(atac)

atac[atac > 1] = 1
colData=colnames(atac)
rowData=rownames(atac)
rowRanges=StringToGRanges(rowData)
rse = SummarizedExperiment::SummarizedExperiment(assays = list(counts = atac),
                                                 rowRanges = rowRanges, colData = colData)

colData(rse)
rowRanges(rse)
assay(rse)[1:10,1:10]



# peak annotation
model = "hg38" # choose from hg19, hg38, mm10
rse = annotateRSE(rse, model)

# quality control
rse = doQC(rse, regionSumCutoff = 5, cellSumCutoffSDs = 3)

# Estimate number of clusters
clusterEst = estimateNClusters(rse, nClustersRange = 2:20)
nClusters = clusterEst$nClustersList$logLikeElbow
plotNClusters(clusterEst)
nClusters = 4

# Perform destin clustering
nCores = 2
clusterResults = destinGrid (rse, nClusters = nClusters, nCores = nCores)
clusterResults$cluster
cluster.results=clusterResults$cluster$cluster

# Visualization by tSNE
PCs = clusterResults$PCs
tsne = Rtsne(as.matrix(PCs))
componentX = tsne$Y[, 1]
componentY = tsne$Y[, 2]

toplot=data.frame(tsne1=componentX, tsne2=componentY,
                  group=hiv$group[names(cluster.results)])
sp<-ggplot(toplot, aes(x=tsne1, y=tsne2, color=group)) + geom_point(size = 0.5)+
  guides(colour = guide_legend(override.aes = list(size=5)))

pdf("pdf/two/figure2.pdf", width = 15, height = 15)
sp
dev.off()


use_python('/opt/anaconda3/bin/python3.9', required = TRUE)

k=2; distance = "euclidean"; n_neighbors = 10; min_dist = 0.1; rand.seed = 42
set.seed(rand.seed)
reticulate::py_set_seed(rand.seed)


UMAP <- reticulate::import("umap")
umapper <- UMAP$UMAP(n_components = as.integer(k), metric = distance, 
                     n_neighbors = as.integer(n_neighbors), min_dist = min_dist)
Rumap <- umapper$fit_transform
umap <- Rumap(as.matrix(PCs))
componentX = umap[, 1]
componentY = umap[, 2]

toplot=data.frame(destinUMAP1=componentX, destinUMAP2=componentY,
                  group=hiv$group[names(cluster.results)])
sp<-ggplot(toplot, aes(x=destinUMAP1, y=destinUMAP2, color=group)) + geom_point(size = 0.5)+
  guides(colour = guide_legend(override.aes = list(size=5)))
sp

ggsave(sp, file='pdf/two/figure3.pdf', width=5, height=5)


DimPlot(hiv, reduction = "umap.atac", group.by = "group", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("scATAC-seq group")

toplot=data.frame(signacUMAP1=hiv@reductions$umap.atac@cell.embeddings[,1],
                  signacUMAP2=hiv@reductions$umap.atac@cell.embeddings[,2],
                  group=hiv$group)
sp<-ggplot(toplot, aes(x=signacUMAP1, y=signacUMAP2, color=group)) + geom_point(size = 0.5)+
  guides(colour = guide_legend(override.aes = list(size=5)))
sp

ggsave(sp, file='pdf/two/figure4.pdf', width=5, height=4)

DefaultAssay(hiv) <- "ATAC"
get_motif_umap=function(TF,...){
  motif.name <- ConvertMotifID(hiv, name = TF)
  motif_plot <- FeaturePlot(hiv, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"))
  motif_plot <- motif_plot + ggtitle(paste('chromVAR access. deviation: TF',TF,'motif', names(which(hiv@assays$ATAC@motifs@motif.names==TF))))
  motif_plot
}
pdf("pdf/two/figure1.pdf", width = 20, height = 5)

get_motif_umap('JUN(var.2)') | get_motif_umap('FOS') | DimPlot(hiv, reduction = "umap.atac",  group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")+ NoLegend()

dev.off()


load("hiv.gene.activity.rda")
pdf("pdf/two/figure5.pdf", width = 10, height = 10)
DimPlot(hiv, reduction = "umap.geneactivity", group.by='group')+ggtitle('Gene Activity')

dev.off()



pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

hiv@assays$ATAC=subset(hiv@assays$ATAC, rownames(hiv)[as.character(seqnames(granges(hiv)))!='HIV'])

# add motif information
hiv<- AddMotifs(hiv, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)
Fragments(hiv@assays$ATAC) <- NULL
fragments <- CreateFragmentObject(path = "/Users/mwen/Desktop/cell line/atac_fragments.tsv.gz", cells = colnames(hiv), validate.fragments = TRUE)
fragments <- CreateFragmentObject(path = "/Users/mwen/Documents/Dissertation/HIV/HIV1/outs/atac_fragments.tsv.gz", cells = colnames(hiv), validate.fragments = TRUE)

Fragments(hiv@assays$ATAC) <- fragments
hiv <- Footprint(
  object = hiv,
  motif.name = c("SNAI2", "SREBF1", "ZNF135", "NFATC2", "RBPJ", "TEAD2", "JUN", "FOS"),
  #motif.name =c("SNAI1", "SNAI2","SNAI3","FIGLA", "JUN", "FOS"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)

hiv <- Footprint(
  object = hiv,
  motif.name = "SNAI2",
  #motif.name =c("SNAI1", "SNAI2","SNAI3","FIGLA", "JUN", "FOS"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)

p2 <- PlotFootprint(hiv, features = c("SNAI2", "SREBF1", "ZNF135", "NFATC2", "RBPJ", "TEAD2", "JUN", "FOS"))
p2 <- PlotFootprint(hiv, features = c("SNAI1", "SNAI2","SNAI3","FIGLA", "JUN", "FOS"))


pdf("footprint.1.pdf", width =15 , height = 10)
p2
dev.off()

cell.line@assays$ATAC=subset(cell.line@assays$ATAC, rownames(cell.line)[as.character(seqnames(granges(cell.line)))!='HIV'])



cellline =  c("SNAI1", "SNAI2","SNAI3","FIGLA", "JUN", "FOS")
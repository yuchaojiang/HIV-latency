setwd("~/Dropbox/HIV")
setwd("C:/Users/yuchaoj/Dropbox/HIV")

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


load('processed_data/hiv_processed.rda')
load('processed_data/transcripts.rda')
DefaultAssay(hiv)='ATAC'

hiv.genes=tail(rownames(hiv@assays$SCT@scale.data))
DotPlot(hiv, features=paste0('sct_', hiv.genes), group.by='group')

hiv.peaks=GRangesToString(hiv@assays$ATAC@ranges[seqnames(hiv@assays$ATAC@ranges)=='HIV'])
library(cowplot)
source('DotPlot.peak.R')
p1=DotPlot.peak(hiv, features=hiv.peaks[1], group.by='group')
p2=DotPlot.peak(hiv, features=hiv.peaks[2], group.by='group')
p3=DotPlot.peak(hiv, features=hiv.peaks[3], group.by='group')
p4=DotPlot.peak(hiv, features=hiv.peaks[4], group.by='group')
p1+p2+p3+p4
ggsave(p3, file='output/supp_figure4_dotplot.peak.pdf', width=4, height=3.5)

p1=DotPlot.peak(hiv, features=hiv.peaks[1], group.by='seurat_clusters')
p2=DotPlot.peak(hiv, features=hiv.peaks[2], group.by='seurat_clusters')
p3=DotPlot.peak(hiv, features=hiv.peaks[3], group.by='seurat_clusters')
p4=DotPlot.peak(hiv, features=hiv.peaks[4], group.by='seurat_clusters')
p1+p2+p3+p4


# Get the total number of HIV ATAC reads
frag.file <- "../Multiome_Shared/data/HIV/outs/atac_fragments.tsv.gz"
DefaultAssay(hiv)='ATAC'
frag.hiv <- CreateFragmentObject(
  path = frag.file,
  cells = colnames(hiv)
)

hiv.ranges= append(hiv@assays$ATAC@ranges[which(as.factor(seqnames(hiv@assays$ATAC@ranges))=='HIV')] # four peaks
                   , GRanges(seqnames='HIV', ranges=IRanges(start=1, end=6003))) # HIV whole chromosome
# The FeatureMatrix function works for at least two peaks
nATAC <- FeatureMatrix(
  fragments = frag.hiv,
  features =hiv.ranges,
  cells = colnames(hiv)
)
all(hiv$hiv_ATAC_total==nATAC[5,])

DotPlot.peak(hiv, features='hiv_ATAC_total', group.by='group', scale = FALSE)
DotPlot.peak(hiv, features='hiv_ATAC_total', group.by='seurat_clusters', scale=FALSE)


hist(hiv$pRNA.sqrt)
par(mfrow=c(2,2))
hist(hiv$pRNA.sqrt[hiv$group=='DMSO'], main='DMSO HIV pRNA', xlab='Proportion of HIV RNA reads', breaks=seq(0,0.2,0.01), xlim=c(0,0.15))
hist(hiv$pRNA.sqrt[hiv$group=='iBET151'], main='iBET151 HIV pRNA', xlab='Proportion of HIV RNA reads', breaks=seq(0,0.2,0.01), xlim=c(0,0.15))
hist(hiv$pRNA.sqrt[hiv$group=='Prostratin'], main='Prostratin HIV pRNA', xlab='Proportion of HIV RNA reads', breaks=seq(0,0.2,0.01), xlim=c(0,0.15))
hist(hiv$pRNA.sqrt[hiv$group=='SAHA'], main='SAHA HIV pRNA', xlab='Proportion of HIV RNA reads', breaks=seq(0,0.2,0.01), xlim=c(0,0.15))

# Test whether HIV exp cells have higher ATAC than HIV non.exp cells
group.temp=as.character(hiv$group)
group.temp[hiv$pRNA.sqrt>0]='HIV.expr'
group.temp[hiv$pRNA.sqrt==0]='HIV.no.expr'
table(group.temp)
group.temp=as.factor(group.temp)


# MEAN
nATAC.df=data.frame(t(nATAC), group=group.temp)

obs=aggregate(nATAC.df[, 1:5], list(nATAC.df$group), mean)
obs
obs=as.numeric(obs[1,2:ncol(obs)]-obs[2,2:ncol(obs)])
obs

nsamp=10000
samp.stat=matrix(nrow=nsamp, ncol=length(obs))
colnames(samp.stat)=c(paste0('nuc-',0:3),'HIV_ATAC_total')
for(i in 1:nsamp){
  if(i %%100 ==0 ) cat(i,'\t')
  group.temp.samp=group.temp[sample(length(group.temp))]
  nATAC.df.samp=data.frame(t(nATAC), group=group.temp.samp)
  samp=aggregate(nATAC.df.samp[, 1:5], list(nATAC.df.samp$group), mean)
  samp.stat[i,]=as.numeric(samp[1,2:ncol(samp)]-samp[2,2:ncol(samp)])
}

pval=rep(NA, length(obs))
for(i in 1:length(pval)){
  pval[i] = mean(samp.stat[,i]>obs[i])
}
par(mfrow=c(2,3))
for(i in 1:5){
  hist(samp.stat[,i], main=paste(colnames(samp.stat)[i], 'p-value =', pval[i]),
       xlab='Test statistic')
  abline(v=obs[i], lty=2, lwd=2, col='red')
}


# PROPORTION of NON-ZERO
nATAC.df=data.frame(t(nATAC), group=group.temp)
obs=aggregate(nATAC.df[, 1:5], list(nATAC.df$group), function(x){sum(x>0)/length(x)})
obs
obs=as.numeric(obs[1,2:ncol(obs)]-obs[2,2:ncol(obs)])
obs

nsamp=10000
samp.stat=matrix(nrow=nsamp, ncol=length(obs))
colnames(samp.stat)=c(paste0('nuc-',0:3),'HIV_ATAC_total')
for(i in 1:nsamp){
  if(i %%100 ==0 ) cat(i,'\t')
  group.temp.samp=group.temp[sample(length(group.temp))]
  nATAC.df.samp=data.frame(t(nATAC), group=group.temp.samp)
  samp=aggregate(nATAC.df.samp[, 1:5], list(nATAC.df.samp$group), function(x){sum(x>0)/length(x)})
  samp.stat[i,]=as.numeric(samp[1,2:ncol(samp)]-samp[2,2:ncol(samp)])
}


pval=rep(NA, length(obs))
for(i in 1:length(pval)){
  pval[i] = mean(samp.stat[,i]>obs[i])
}
par(mfrow=c(2,3))
for(i in 1:5){
  hist(samp.stat[,i], main=paste(colnames(samp.stat)[i], 'p-value =', pval[i]),
       xlab='Test statistic')
  abline(v=obs[i], lty=2, lwd=2, col='red')
}


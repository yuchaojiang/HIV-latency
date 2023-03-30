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


load('processed_data/hiv.rda')
load('processed_data/transcripts.rda')
DefaultAssay(hiv)='RNA'

# Treatment group
table(hiv$group)

# percentage of RNA counts from HIV 
pRNA=  apply(hiv@assays$RNA@counts[which(as.factor(seqnames(transcripts))=='HIV'),],2,sum)/
  apply(hiv@assays$RNA@counts,2,sum)
boxplot(sqrt(pRNA)~hiv$group, xlab='', ylab='sqrt(% of RNA reads from HIV)')  
  
hiv$pRNA.sqrt=sqrt(pRNA)
VlnPlot(hiv, features='pRNA.sqrt')+ggtitle('sqrt(% of RNA reads from HIV)')


# percentage of ATAC counts from HIV
pATAC= apply(hiv@assays$ATAC@counts[which(as.factor(seqnames(hiv@assays$ATAC@ranges))=='HIV'),],2,sum)/
  apply(hiv@assays$ATAC@counts,2,sum)
boxplot(sqrt(pATAC)~hiv$group)# Does not seem to make sense
plot(sqrt(pRNA), sqrt(pATAC), col=hiv$group)

# Total number of HIV ATAC reads within the peaks: at most 3 ATAC reads
# Total counts from the peak matrix
table(apply(hiv@assays$ATAC@counts[which(as.factor(seqnames(hiv@assays$ATAC@ranges))=='HIV'),],2,sum))
# Reprocess
frag.file <- "../Multiome_Shared/data/HIV/outs/atac_fragments.tsv.gz"
DefaultAssay(hiv)='ATAC'
frag.hiv <- CreateFragmentObject(
  path = frag.file,
  cells = colnames(hiv)
)

hiv.ranges= append(hiv@assays$ATAC@ranges[which(as.factor(seqnames(hiv@assays$ATAC@ranges))=='HIV')] # four peaks
                   , GRanges(seqnames='HIV', ranges=IRanges(start=1, end=6003))) # HIV whole chromosome
# The FeatureMatrix function works for at least two peaks
nATAC_total <- FeatureMatrix(
  fragments = frag.hiv,
  features =hiv.ranges,
  cells = colnames(hiv)
)
hiv_ATAC_total=nATAC_total[5,]
boxplot( hiv_ATAC_total~ hiv$group, xlab='', ylab='Total HIV ATAC reads')

hist(hiv_ATAC_total,
     breaks=seq(-0.5, 3.5,1), xlab='Total reads from HIV', ylab='No. of cells',
     main='Total number of HIV ATAC reads')
hiv$hiv_ATAC_total=hiv_ATAC_total

# third HIV ATAC peak (nuc-3)
ATAC_count_nuc3=hiv@assays$ATAC@counts[nrow(hiv@assays$ATAC@counts)-1,]
hist(ATAC_count_nuc3,breaks=seq(-0.5, 5.5,1))
table(ATAC_count_nuc3)

barplot(table(ATAC_count_nuc3,hiv$group)[2,]/table(hiv$group),
        ylab='Percentage of Open-Chromatin Cells in Nuc-3')

smoothScatter(hiv_ATAC_total, sqrt(pRNA), cex=1, pch=16, nrpoints = 0,
              xlab='Total HIV ATAC reads', ylab='sqrt percentage of HIV RNA reads')
points(hiv_ATAC_total, sqrt(pRNA), pch=16, cex=0.5, col='gray')
lm.run=lm(sqrt(pRNA)~hiv_ATAC_total)
abline(lm.run, col='black', lty=2)
title(paste('coeff =', signif(summary(lm.run)$coefficients[2,1],3), 
            'pval =', signif(summary(lm.run)$coefficients[2,4],3)))


par(mfrow=c(2,2))
smoothScatter(nATAC_total[1,], sqrt(pRNA), cex=1, pch=16, nrpoints = 0,
              xlab='Nuc-0 ATAC reads', ylab='sqrt percentage of HIV RNA reads')
points(nATAC_total[1,], sqrt(pRNA), pch=16, cex=0.5, col='gray')
lm.run=lm(sqrt(pRNA)~nATAC_total[1,])
abline(lm.run, col='black', lty=2)
title(paste('coeff =', signif(summary(lm.run)$coefficients[2,1],3), 
            'pval =', signif(summary(lm.run)$coefficients[2,4],3)))

smoothScatter(nATAC_total[2,], sqrt(pRNA), cex=1, pch=16, nrpoints = 0,
              xlab='Nuc-1 ATAC reads', ylab='sqrt percentage of HIV RNA reads')
points(nATAC_total[2,], sqrt(pRNA), pch=16, cex=0.5, col='gray')
lm.run=lm(sqrt(pRNA)~nATAC_total[2,])
abline(lm.run, col='black', lty=2)
title(paste('coeff =', signif(summary(lm.run)$coefficients[2,1],3), 
            'pval =', signif(summary(lm.run)$coefficients[2,4],3)))

table(nATAC_total[3,])
smoothScatter(nATAC_total[3,], sqrt(pRNA), cex=1, pch=16, nrpoints = 0,
              xlab='Nuc-2 ATAC reads', ylab='sqrt percentage of HIV RNA reads')
points(nATAC_total[3,], sqrt(pRNA), pch=16, cex=0.5, col='gray')
lm.run=lm(sqrt(pRNA)~nATAC_total[3,])
abline(lm.run, col='black', lty=2)
title(paste('coeff =', signif(summary(lm.run)$coefficients[2,1],3), 
            'pval =', signif(summary(lm.run)$coefficients[2,4],3)))

smoothScatter(nATAC_total[4,], sqrt(pRNA), cex=1, pch=16, nrpoints = 0,
              xlab='3\' peak ATAC reads', ylab='sqrt percentage of HIV RNA reads')
points(nATAC_total[4,], sqrt(pRNA), pch=16, cex=0.5, col='gray')
lm.run=lm(sqrt(pRNA)~nATAC_total[4,])
abline(lm.run, col='black', lty=2)
title(paste('coeff =', signif(summary(lm.run)$coefficients[2,1],3), 
            'pval =', signif(summary(lm.run)$coefficients[2,4],3)))



par(mfrow=c(1,1))
save(hiv, file='processed_data/hiv_processed.rda')

# Below is to generate input for Alex and Cynthia
# save(transcripts, 'to_share/transcripts.rda')
# 
# atac.mat=hiv@assays$ATAC@counts
# save(atac.mat, file='to_share/atac.mat.rda')
# 
# atac.ranges=hiv@assays$ATAC@ranges
# save(atac.ranges, file='to_share/atac.ranges.rda')
# 
# rna.mat=hiv@assays$RNA@counts
# save(rna.mat, file='to_share/rna.mat.rda')
# rm(rna.mat)
# 
# rna.normalized.mat=hiv@assays$SCT@data
# save(rna.normalized.mat, file='to_share/rna.normalized.mat.rda')
# rm(rna.normalized.mat)
# 
# group=hiv$group
# save(group, file='to_share/group.rda')
# rm(group)

# motif score

load('processed_data/hiv_chromvar.rda')
motif.mat=hiv@assays$chromvar@data
TF=ConvertMotifID(hiv, id=rownames(motif.mat))
# save(motif.mat, file='to_share/motif.mat.rda')
# save(TF, file='to_share/TF.rda')

load('processed_data/hiv_processed.rda')

par(mfrow=c(1,3))
TF.name='JUN(var.2)'
#boxplot(motif.mat[which(TF==TF.name),]~hiv$group)
temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
smoothScatter(motif.mat[which(TF==TF.name),], sqrt(pRNA), xlab='TF motif score',
              ylab=paste('sqrt(percentage of HIV RNA reads)'),
              main=paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                         ': p-val =',signif(temp$coefficients[2,4],4)))
abline(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]), col='black',lty=2)


TF.name='FOS'
#boxplot(motif.mat[which(TF==TF.name),]~hiv$group)
temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
smoothScatter(motif.mat[which(TF==TF.name),], sqrt(pRNA), xlab='TF motif score',
              ylab=paste('sqrt(percentage of HIV RNA reads)'),
              main=paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                         ': p-val =',signif(temp$coefficients[2,4],4)))
abline(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]), col='black',lty=2)

TF.name='FOS::JUN'
#boxplot(motif.mat[which(TF==TF.name),]~hiv$group)
temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
smoothScatter(motif.mat[which(TF==TF.name),], sqrt(pRNA), xlab='TF motif score',
              ylab=paste('sqrt(percentage of HIV RNA reads)'),
              main=paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                         ': p-val =',signif(temp$coefficients[2,4],4)))
abline(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]), col='black',lty=2)
par(mfrow=c(1,1))

# Get the p-value for each TF
linked.TF=matrix(nrow=length(TF), ncol=4)
linked.TF=as.data.frame(linked.TF)
linked.TF[,1]=TF
linked.TF[,2]=rownames(motif.mat)
colnames(linked.TF)=c('TF','motif','cor','pval')
for(i in 1:length(TF)){
  TF.name=TF[i]
  temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
  linked.TF[i,3]=signif(cor(sqrt(pRNA), motif.mat[which(TF==TF.name),], method='spearman'),4)
  linked.TF[i,4]=signif(temp$coefficients[2,4],4)
}
pdf(file='output/figure5_pval_link_TF.pdf', width=4, height=4)
hist(linked.TF[,4], xlab='p-value', main='P-val dist. from TF linkage analysis',
     breaks=seq(0,1,0.05))
dev.off()

library(qvalue)
linked.TF=cbind(linked.TF, pval.adj=signif(qvalue(p=linked.TF[,4])$qvalues,4))

pdf(file='output/figure5_pval_link_TF_volcano.pdf', width=4, height=4)
plot(linked.TF[,3], -log(linked.TF[,5],10), pch=16, xlab='Spearman correlation', ylab='-log10(adj. p-value)', cex=0.5)
sel=(-log(linked.TF[,5],10))>=(-log(0.05,10))
points(linked.TF[sel,3], -log(linked.TF[sel,5],10), pch=16, col=2, cex=0.5)
abline(h=-log(0.05,10), lty=2)
title('TF linkage analysis')
dev.off()

head(linked.TF[order(linked.TF[,5]),])


write.csv(linked.TF[order(linked.TF[,5]),], file='linkage_analysis_output/linked.TF.csv', row.names = F)
temp=linked.TF[order(linked.TF[,5]),]
head(temp)
temp=temp[temp$pval.adj<=0.05,]
table(temp$cor>0)
quantile(temp$cor^2, c(0.025, 0.975))
paste(temp$TF[temp$cor>0], collapse = ' ')
paste(temp$TF[temp$cor<0], collapse = ' ')



# Plot for top linked TFs
pdf(file='output/figure5_pval_link_TF_top_scatter.pdf', width=9, height=3)
par(mfrow=c(1,3))
for(TF.name in c('SNAI2','ELF1','GABPA')){
  temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
  smoothScatter(motif.mat[which(TF==TF.name),],sqrt(pRNA), xlab='TF motif score.', ylab='sqrt(percentage of HIV RNA reads)',
                main=paste(TF.name, '\nr =',signif(cor(sqrt(pRNA), motif.mat[which(TF==TF.name),], method='spearman'),3),
                           'p-val =',signif(temp$coefficients[2,4],3)))
  abline(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]), lty=2, col='black')
  
  # # below is applying regression to only cells with pRNA>0
  # temp=summary(lm(sqrt(pRNA)[pRNA>0]~motif.mat[which(TF==TF.name),][pRNA>0]))
  # smoothScatter(motif.mat[which(TF==TF.name),][pRNA>0],sqrt(pRNA)[pRNA>0], xlab='TF motif score.', ylab='sqrt(percentage of HIV RNA reads)',
  #               main=paste(TF.name, 'p-val =',signif(temp$coefficients[2,4],4)))
  # abline(lm(sqrt(pRNA)[pRNA>0]~motif.mat[which(TF==TF.name),][pRNA>0]), lty=2, col='black')
  
}
par(mfrow=c(1,1))
dev.off()



# Get the p-value for each gene
var.genes=hiv@assays$SCT@var.features

linked.gene=matrix(nrow=length(var.genes), ncol=3)
linked.gene=as.data.frame(linked.gene)
linked.gene[,1]=var.genes
colnames(linked.gene)=c('gene','cor','pval')
for(i in 1:nrow(linked.gene)){
  gene.name=var.genes[i]
  if(i %%200==0) cat(i,'\t')
  temp=summary(lm(sqrt(pRNA)~hiv@assays$SCT@data[gene.name,]))
  linked.gene[i,2]=signif(cor(sqrt(pRNA), hiv@assays$SCT@data[gene.name,], method='spearman'),4)
  linked.gene[i,3]=signif(temp$coefficients[2,4],4)
}

pdf(file='output/figure5_pval_link_gene.pdf', width=4, height=4)
hist(linked.gene[,3], xlab='p-value', main='P-val dist. from gene linkage analysis',
     breaks=seq(0,1,0.05))
dev.off()

library(qvalue)
linked.gene=cbind(linked.gene, pval.adj=signif(qvalue(p=linked.gene[,3])$qvalues,4))

pdf(file='output/figure5_pval_link_gene_volcano.pdf', width=4, height=4)
plot(linked.gene[,2], -log(linked.gene[,4],10), pch=16, xlab='Spearman correlation', ylab='-log(adj. p-value)', cex=0.5,
     ylim=c(0,20), xlim=c(-0.15, 0.15))
sel=(-log(linked.gene[,4],10))>=(-log(0.05,10))
points(linked.gene[sel,2], -log(linked.gene[sel,4],10), pch=16, col=2, cex=0.5)
abline(h=-log(0.05,10), lty=2)
title('Gene linkage analysis')
dev.off()

linked.gene.output=linked.gene[order(linked.gene[,4]),]
write.csv(linked.gene.output, file='linkage_analysis_output/linked.gene.csv', row.names = F)
linked.gene.output[1:20,]

# Plot for top linked genes
pdf(file='output/figure5_pval_link_gene_top_scatter.pdf', width=9, height=3)
par(mfrow=c(1,3))
for(gene.name in c('TSPOAP1','MALAT1','NEAT1')){
  temp=summary(lm(sqrt(pRNA)~hiv@assays$SCT@data[gene.name,]))
  smoothScatter(hiv@assays$SCT@data[gene.name,],sqrt(pRNA), xlab='RNA exp.', ylab='sqrt(percentage of HIV RNA reads)',
                main=paste(gene.name, '\nr =', signif(cor(sqrt(pRNA), hiv@assays$SCT@data[gene.name,], method='spearman'),3),
                           'p-val =',signif(temp$coefficients[2,4],3)))
  abline(lm(sqrt(pRNA)~hiv@assays$SCT@data[gene.name,]), lty=2, col='black')
}
par(mfrow=c(1,1))
dev.off()


# Get the p-value for each TF gene expression
TF.genes=TF[TF %in% rownames(hiv@assays$SCT)]

linked.TF.gene=matrix(nrow=length(TF.genes), ncol=3)
linked.TF.gene=as.data.frame(linked.TF.gene)
linked.TF.gene[,1]=TF.genes
colnames(linked.TF.gene)=c('gene','cor','pval')
for(i in 1:nrow(linked.TF.gene)){
  gene.name=TF.genes[i]
  if(i %%200==0) cat(i,'\t')
  temp=summary(lm(sqrt(pRNA)~hiv@assays$SCT@data[gene.name,]))
  linked.TF.gene[i,2]=signif(cor(sqrt(pRNA), hiv@assays$SCT@data[gene.name,], method='spearman'),4)
  linked.TF.gene[i,3]=signif(temp$coefficients[2,4],4)
}

pdf(file='output/figure5_pval_link_TF_exp.pdf', width=4, height=4)
hist(linked.TF.gene[,3], xlab='p-value', main='P-val dist. from TF exp. linkage analysis',
     breaks=seq(0,1,0.05))
dev.off()

library(qvalue)
linked.TF.gene=cbind(linked.TF.gene, pval.adj=signif(qvalue(p=linked.TF.gene[,3])$qvalues,4))

pdf(file='output/figure5_pval_link_TF_exp_volcano.pdf', width=4, height=4)
plot(linked.TF.gene[,2], -log(linked.TF.gene[,4],10), pch=16, xlab='Spearman correlation', ylab='-log(adj. p-value)', cex=0.5,
     ylim=c(0,4), xlim=c(-0.07, 0.07))
sel=(-log(linked.TF.gene[,4],10))>=(-log(0.05,10))
points(linked.TF.gene[sel,2], -log(linked.TF.gene[sel,4],10), pch=16, col=2, cex=0.5)
abline(h=-log(0.05,10), lty=2)
title('TF exp. linkage analysis')
dev.off()

linked.TF.gene.output=linked.TF.gene[order(linked.TF.gene[,4]),]
colnames(linked.TF.gene.output)[1]='TF.gene'
head(linked.TF.gene.output)
write.csv(linked.TF.gene.output, file='linkage_analysis_output/linked.TF.gene.csv', row.names = F)

# Plot for top linked TF genes
pdf(file='output/figure5_pval_link_TF_exp_top_scatter.pdf', width=9, height=3)
par(mfrow=c(1,3))
for(gene.name in c('GATA3','FOXA3','ELF1')){
  temp=summary(lm(sqrt(pRNA)~hiv@assays$SCT@data[gene.name,]))
  smoothScatter(hiv@assays$SCT@data[gene.name,],sqrt(pRNA), xlab='TF RNA exp.', ylab='sqrt(percentage of HIV RNA reads)',
                main=paste(gene.name, '\nr =', signif(cor(sqrt(pRNA), hiv@assays$SCT@data[gene.name,], method='spearman'),3),
                           'p-val =',signif(temp$coefficients[2,4],3)))
  abline(lm(sqrt(pRNA)~hiv@assays$SCT@data[gene.name,]), lty=2, col='black')
}
par(mfrow=c(1,1))
dev.off()


# Get the p-value for each peak
var.peaks=hiv@assays$ATAC@var.features[1:20000]

linked.peak=matrix(nrow=length(var.peaks), ncol=3)
linked.peak=as.data.frame(linked.peak)
linked.peak[,1]=var.peaks
colnames(linked.peak)=c('peak','cor','pval')
for(i in 1:nrow(linked.peak)){
  peak.index=which(rownames(hiv@assays$ATAC)==var.peaks[i])
  if(i %%400==0) cat(i,'\t')
  temp=summary(lm(sqrt(pRNA)~hiv@assays$ATAC@data[peak.index,]))
  linked.peak[i,2]=signif(cor(sqrt(pRNA), hiv@assays$ATAC@data[peak.index,], method='spearman'),4)
  linked.peak[i,3]=signif(temp$coefficients[2,4],4)
}


pdf(file='output/figure5_pval_link_peak.pdf', width=4, height=4)
hist(linked.peak[,3], xlab='p-value', main='P-val dist. from peak linkage analysis',
     breaks=seq(0,1,0.05))
dev.off()

# # There are too many peaks, and FDR is too stringent.
# library(qvalue)
# linked.peak=cbind(linked.peak, pval.adj=qvalue(p=linked.peak[,3])$qvalues)

#  We will just use nominal p-values

pdf(file='output/figure5_pval_link_peak_volcano.pdf', width=4, height=4)
plot(linked.peak[,2], -log(linked.peak[,3],10), pch=16, xlab='Spearman correlation', ylab='-log(p-value)', cex=0.5,
     ylim=c(0,4), xlim=c(-0.1, 0.1))
sel=(-log(linked.peak[,3],10))>=(-log(0.01,10))
points(linked.peak[sel,2], -log(linked.peak[sel,3],10), pch=16, col=2, cex=0.5)
abline(h=-log(0.01,10), lty=2)
title('Peak linkage analysis')
dev.off()



linked.peak.output=linked.peak[order(linked.peak[,3]),] # only output significant ones
write.csv(linked.peak.output, file='linkage_analysis_output/linked.peak.csv', row.names = F)

load('processed_data/hiv_chromvar.rda')
seqinfo(hiv)=seqinfo(hiv)[seqnames(hiv)[1:25]]

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
hiv <- AddMotifs(
  object = hiv,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# All peaks
sig.link.peaks= linked.peak.output[linked.peak.output[,3]<=0.01,1]
sig.link.peaks=sig.link.peaks[!grepl('HIV', sig.link.peaks)] 
sig.link.peaks=sig.link.peaks[!grepl('MT', sig.link.peaks)]

set.seed(1)
enriched.motifs <- FindMotifs(
  object = hiv,
  features =sig.link.peaks
)
enriched.motifs
write.csv(enriched.motifs, file='linkage_analysis_output/enriched.motif.csv', row.names = F)
p=MotifPlot(
  object = hiv,
  motifs = head(rownames(enriched.motifs)),
  ncol = 6
)
p

# Positively linked peaks
sig.pos.link.peaks= linked.peak.output[linked.peak.output$pval<=0.01 & linked.peak.output$cor>0,1]
sig.pos.link.peaks=sig.pos.link.peaks[!grepl('HIV', sig.pos.link.peaks)] 
sig.pos.link.peaks=sig.pos.link.peaks[!grepl('MT', sig.pos.link.peaks)]

set.seed(1)
pos.enriched.motifs <- FindMotifs(
  object = hiv,
  features =sig.pos.link.peaks
)
pos.enriched.motifs
write.csv(pos.enriched.motifs, file='linkage_analysis_output/pos.enriched.motif.csv', row.names = F)
p=MotifPlot(
  object = hiv,
  motifs = head(rownames(pos.enriched.motifs)),
  ncol = 6
)
p
ggsave(p, file='output/figure5_pos_enriched_motif.pdf', width=11, height=1.5)

# Negatively linked peaks
sig.neg.link.peaks= linked.peak.output[linked.peak.output$pval<=0.01 & linked.peak.output$cor<0,1]
sig.neg.link.peaks=sig.neg.link.peaks[!grepl('HIV', sig.neg.link.peaks)] 
sig.neg.link.peaks=sig.neg.link.peaks[!grepl('MT', sig.neg.link.peaks)]

set.seed(1)
neg.enriched.motifs <- FindMotifs(
  object = hiv,
  features =sig.neg.link.peaks
)
neg.enriched.motifs
write.csv(neg.enriched.motifs, file='linkage_analysis_output/neg.enriched.motif.csv', row.names = F)
p=MotifPlot(
  object = hiv,
  motifs = head(rownames(neg.enriched.motifs)),
  ncol = 6
)
p
ggsave(p, file='output/figure5_neg_enriched_motif.pdf', width=11, height=1.5)


library(ggvenn)
x=list(linked.TF=linked.TF[linked.TF$pval.adj<=0.1,2],
       enriched.motif=enriched.motifs[enriched.motifs$pvalue<=0.01,1])
p1=ggvenn(x, 
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)
p1

overlapped.motifs=intersect(linked.TF[linked.TF$pval.adj<=0.1,2], enriched.motifs[enriched.motifs$pvalue<=0.01,1])
linked.TF[match(overlapped.motifs, linked.TF[,2]),]
enriched.motifs[match(overlapped.motifs, enriched.motifs[,1]),]

load('processed_data/hiv.rda')
hiv$pRNA=pRNA
hiv$pRNA.sqrt=sqrt(pRNA)
hiv$hiv_ATAC_total=hiv_ATAC_total

p1 <- DimPlot(hiv, reduction = "umap.rna", group.by='group', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") + NoLegend()
p2 <- DimPlot(hiv, reduction = "umap.atac", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")+ NoLegend()
p3 <- DimPlot(hiv, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p4 <- DimPlot(hiv, reduction = "wnn.umap", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3 +p4 +plot_layout(ncol=2)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))


p1 <- DimPlot(hiv, reduction = "umap.rna", group.by='group', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") + NoLegend()
p2 <- DimPlot(hiv, reduction = "umap.atac", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")+ NoLegend()
p3 <- FeaturePlot(hiv, features = 'pRNA.sqrt', reduction = "umap.rna")
p4 <- FeaturePlot(hiv, features = 'pRNA.sqrt', reduction = "umap.atac")
p1 + p2 + p3 +p4 +plot_layout(ncol=2)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))


p1 <- DimPlot(hiv, reduction = "umap.rna", group.by='group', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") + NoLegend()
p2 <- DimPlot(hiv, reduction = "umap.rna", group.by='seurat_clusters', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") + NoLegend()
p3<-VlnPlot(hiv, features = 'pRNA.sqrt')
p4<-VlnPlot(hiv, features = 'pRNA.sqrt', group.by = 'seurat_clusters')

p1 + p2 + p3 +p4 +plot_layout(ncol=2)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))

cov_plot1 <- CoveragePlot(
  object = hiv,
  region = "HIV-1-6000",
  annotation = FALSE,
  peaks = TRUE
)
cov_plot2 <- CoveragePlot(
  object = hiv,
  region = "HIV-1-6000",
  annotation = FALSE,
  peaks = TRUE,
  group.by='seurat_clusters'
)
cov_plot1|cov_plot2

save.image(file='processed_data/hiv_predict.rda')


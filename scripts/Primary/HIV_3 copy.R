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
library(car)

load("hiv_chromvar.rda")
hiv 
# 48673 features across 61708 samples within 4 assays 
motif.mat=hiv@assays$chromvar@data
save(motif.mat, file='motif.mat.rda')

TF=ConvertMotifID(hiv,  id=rownames(motif.mat))
save(TF, file='TF.rda')

load('hiv.rda')
load('transcripts.rda')
DefaultAssay(hiv)='RNA'
pRNA <- hiv@assays$RNA@counts[which(as.factor(seqnames(transcripts))=='HIV'),]/apply(hiv@assays$RNA@counts,2,sum)
save(pRNA, file="pRNA.rda")
pdf("pdf/three/figure1.pdf", width = 10, height = 10)
boxplot(sqrt(pRNA)~hiv$group, xlab='', ylab='sqrt(% of RNA reads from HIV)')
dev.off()

hiv$pRNA.sqrt=sqrt(pRNA)
hiv$pRNA=pRNA
load("pRNA.rda")


frag.file <- "new.fragment1.txt.gz"
DefaultAssay(hiv)='ATAC'
frag.hiv <- CreateFragmentObject(
  path = frag.file,
  cells = colnames(hiv)
)

hiv.ranges= append(hiv@assays$ATAC@ranges[which(as.factor(seqnames(hiv@assays$ATAC@ranges))=='HIV')], GRanges(seqnames='HIV', ranges=IRanges(start=1, end=10425))) # HIV whole chromosome

# percentage of ATAC counts from HIV
p1= apply(hiv@assays$ATAC@counts[which(as.factor(seqnames(hiv@assays$ATAC@ranges))=='HIV'),],2,sum)
p2 = apply(hiv@assays$ATAC@counts,2,sum)

pATAC = p1/p2
hiv$pATAC.sqrt = sqrt(pATAC)
save(pATAC, file="pATAC.rda")

pdf("pdf/three/figure2.pdf", width = 5, height = 5)
VlnPlot(hiv, features='pRNA.sqrt', pt.size=0,group.by = "group")+ggtitle('sqrt(% of RNA reads from HIV)')
VlnPlot(hiv, features = 'pATAC.sqrt', pt.size=0,group.by = "group")+ggtitle('sqrt(% of ATAC reads from HIV)')

dev.off()

hiv.genes=transcripts[seqnames(transcripts)=='HIV']$gene_name
pdf("pdf/one/figure3.pdf", width = 5, height = 5)
DotPlot(hiv, features=paste0('sct_', hiv.genes), group.by='group')
#VlnPlot(hiv, features=paste0('sct_', hiv.genes), group.by='group')
dev.off()

temp = as.data.frame(cbind(hiv$pATAC.sqrt,hiv$group))
names(temp) <- c("reads","Group")
temp$Group = as.character(temp$Group)
pdf("pdf/three/figure2.1.pdf", width = 10, height = 10)

ggplot() + 
  geom_boxplot(data = temp, mapping = aes(Group, reads)) + ggtitle("sqrt(% of ATAC reads from HIV)") +
  theme( axis.text.x = element_blank())

dev.off()

anova <- aov(hiv$pATAC.sqrt ~ hiv$group)
summary(anova)
tukey<-TukeyHSD(anova)
tukey

anova <- aov(hiv$pRNA.sqrt ~ hiv$group)
summary(anova)
tukey<-TukeyHSD(anova)
tukey


nATAC_total <- FeatureMatrix(
  fragments = frag.hiv,
  features =hiv.ranges,
  cells = colnames(hiv)
)

hiv_ATAC_total=nATAC_total[3,]
read.count = data.frame(cbind(cbind(pRNA,pATAC),hiv_ATAC_total))

libsize <- colSums(hiv@assays$ATAC@counts)
hiv$accessibility = hiv_ATAC_total/libsize

temp = subset(x = hiv, subset = accessibility >0)
bar = as.data.frame(table(temp$group)/table(hiv$group))
col= c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
pdf("pdf/three/figureadd.pdf",width=10,height = 10)
VlnPlot(hiv, features='accessibility', pt.size=0)+ggtitle('Fraction of Cells that are Open for HIV Chromosome')
barplot(bar$Freq,names.arg=bar$Var1, xlab = "Group", ylab="Fraction of Cells that are Open for HIV Chromosome", col=col)
dev.off()



pdf("pdf/three/figure2.2.pdf", width = 5, height = 5)
VlnPlot(hiv, features='pRNA.sqrt', pt.size=0,group.by = "seurat_clusters")+ggtitle('sqrt(% of RNA reads from HIV)') + theme(legend.position = 'none')
VlnPlot(hiv, features = 'pATAC.sqrt', pt.size=0,group.by = "seurat_clusters")+ggtitle('sqrt(% of ATAC reads from HIV)') + theme(legend.position = 'none')

dev.off()




pdf("pdf/three/figureadd.pdf",width=20,height = 7)
p1 = FeaturePlot(hiv, features = 'accessibility', reduction = "umap.atac",split.by="group",keep.scale="all") & theme(legend.position = "right")
p1
p2 = FeaturePlot(hiv, features = 'accessibility', reduction = "umap.atac",split.by="group",keep.scale="all",cells=colnames(hiv)[hiv$accessibility>0]) & theme(legend.position = "right")
p2
dev.off()

barplot(H,names.arg=M,xlab="Month",ylab="Revenue",col="blue",
        main="Revenue chart",border="red")
  
library("ggpubr")
hiv$pRNA = pRNA
ggscatter(read.count, x = "pRNA", y = "hiv_ATAC_total", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson")
ggscatter(read.count, x = "pRNA", y = "pATAC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson")


smoothScatter(hiv_ATAC_total,sqrt(pRNA),  cex=1, pch=16, nrpoints = 0,
              xlab='Total HIV ATAC reads', ylab='sqrt percentage of HIV RNA reads')
points(hiv_ATAC_total, sqrt(pRNA), pch=16, cex=0.5, col='gray')
lm.run=lm(sqrt(pRNA)~hiv_ATAC_total)
abline(lm.run, col='black', lty=2)
title(paste('coeff =', signif(summary(lm.run)$coefficients[2,1],3), 
            'pval =', signif(summary(lm.run)$coefficients[2,4],3)))

smoothScatter(sqrt(pATAC), sqrt(pRNA), cex=1, pch=16, nrpoints = 0,
              xlab='sqrt percentage of HIV ATAC reads', ylab='sqrt percentage of HIV RNA reads')
points(sqrt(pATAC), sqrt(pRNA), pch=16, cex=0.5, col='gray')
lm.run=lm(sqrt(pRNA)~sqrt(pATAC))
abline(lm.run, col='black', lty=2)
title(paste('coeff =', signif(summary(lm.run)$coefficients[2,1],3), 
            'pval =', signif(summary(lm.run)$coefficients[2,4],3)))




pdf("pdf/three/figure3.pdf", width = 10, height = 10)

boxplot(hiv$pATAC.sqrt~hiv$group, col=c("#F8766D","#7CAE00","#00BFC4","#C77CFF"))
# Does not seem to make sense    
plot(hiv$pRNA.sqrt, hiv$pATAC.sqrt, col=hiv$group)
dev.off()


hiv$hiv_ATAC_total=hiv_ATAC_total

pdf("pdf/three/figureed.pdf", width = 10, height = 10)

p1 <- DimPlot(hiv, reduction = "umap.rna", group.by='group', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv RNA") + NoLegend()
p2 <- DimPlot(hiv, reduction = "umap.atac", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X hiv ATAC")+ NoLegend()
p3 <- FeaturePlot(hiv, features = 'pRNA.sqrt', reduction = "umap.rna")
p4 <- FeaturePlot(hiv, features = 'pATAC.sqrt', reduction = "umap.atac")
p1 + p2 + p3 +p4 +plot_layout(ncol=2)& theme(plot.title = element_text(hjust = 0.5))

dev.off()


pdf("pdf/three/figure4.pdf", width = 10, height = 10)
boxplot(sqrt(pATAC)~ hiv$group, xlab='', ylab='Proportion of HIV ATAC reads')


hist(sqrt(pATAC), xlab='Proportions of HIV ATAC reads', ylab='No. of cells',
     main='Total number of HIV ATAC reads', breaks = 17)


smoothScatter(sqrt(pATAC), sqrt(pRNA), cex=1, pch=16, nrpoints = 0,
              xlab='sqrt percentage of HIV ATAC reads', ylab='sqrt percentage of HIV RNA reads')
points(sqrt(pATAC), sqrt(pRNA), pch=16, cex=0.5, col='gray')
lm.run=lm(sqrt(pRNA)~sqrt(pATAC))
abline(lm.run, col='black', lty=2)
title(paste('coeff =', signif(summary(lm.run)$coefficients[2,1],3),
            'pval =', signif(summary(lm.run)$coefficients[2,4],3)))

dev.off()

save(hiv, file='hiv_processed.rda')


pdf("pdf/three/figure5.pdf", width = 15, height = 5)
color <- c(SAHA="#C77CFF", Prostratin="#00BFC4", iBET151="#7CAE00", DMSO="#F8766D")
colnum <- match(hiv$group, names(color))


par(mfrow=c(1,3))
TF.name='JUN(var.2)'
#boxplot(motif.mat[which(TF==TF.name),]~hiv$group)
temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
smoothScatter(motif.mat[which(TF==TF.name),], sqrt(pRNA),col = color[colnum],  xlab='TF motif score',
              pch = 20, cex = 1, ylab=paste('sqrt(percentage of HIV RNA reads)'),
              main=paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                         ': p-val =',signif(temp$coefficients[2,4],4)))
abline(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]), col='black',lty=2)


TF.name='FOS'
#boxplot(motif.mat[which(TF==TF.name),]~hiv$group)
temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
smoothScatter(motif.mat[which(TF==TF.name),], sqrt(pRNA),col = color[colnum],  xlab='TF motif score',
              pch = 20, cex = 1,  ylab=paste('sqrt(percentage of HIV RNA reads)'),
              main=paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                         ': p-val =',signif(temp$coefficients[2,4],4)))
abline(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]), col='black',lty=2)

TF.name='FOS::JUN'
#boxplot(motif.mat[which(TF==TF.name),]~hiv$group)
temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
smoothScatter(motif.mat[which(TF==TF.name),], sqrt(pRNA),col = color[colnum],  xlab='TF motif score',
              pch = 20, cex = 1, ylab=paste('sqrt(percentage of HIV RNA reads)'),
              main=paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                         ': p-val =',signif(temp$coefficients[2,4],4)))
abline(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]), col='black',lty=2)
par(mfrow=c(1,1))

dev.off()


library(ggplot2)
pdf("pdf/three/figure6.pdf", width = 17, height = 5)

par(mfrow=c(1,3))
TF.name='JUN(var.2)'
#boxplot(motif.mat[which(TF==TF.name),]~hiv$group)
temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
temp1 = data.frame( cbind(motif.mat[which(TF==TF.name),],sqrt(pRNA)))
temp2 = paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
              ': p-val =',signif(temp$coefficients[2,4],4))

p1 = ggplot(temp1, aes(motif.mat[which(TF==TF.name),], sqrt(pRNA))) +
  geom_point(aes(color = hiv$group)) +
  geom_smooth(method='lm') +
  scale_color_manual(values = color) +
  labs(x= 'TF motif score', y = 'sqrt(percentage of HIV RNA reads)', title = paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                                                                                   ': p-val =',signif(temp$coefficients[2,4],4))) +
  theme_bw()+ theme(legend.position = "none")

TF.name='FOS'
#boxplot(motif.mat[which(TF==TF.name),]~hiv$group)
temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
temp1 = data.frame( cbind(motif.mat[which(TF==TF.name),],sqrt(pRNA)))
p2 = ggplot(temp1, aes(motif.mat[which(TF==TF.name),], sqrt(pRNA))) +
  geom_point(aes(color = hiv$group)) +
  geom_smooth(method='lm') +
  scale_color_manual(values = color) +
  labs(x= 'TF motif score', y = 'sqrt(percentage of HIV RNA reads)', title = paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                                                                                   ': p-val =',signif(temp$coefficients[2,4],4))) +
  theme_bw() + theme(legend.position = "none")


TF.name='FOS::JUN'
temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
temp1 = data.frame( cbind(motif.mat[which(TF==TF.name),],sqrt(pRNA)))
p3 = ggplot(temp1, aes(motif.mat[which(TF==TF.name),], sqrt(pRNA))) +
  geom_point(aes(color = hiv$group)) +
  geom_smooth(method='lm') +
  scale_color_manual(values = color) +
  labs(x= 'TF motif score', y = 'sqrt(percentage of HIV RNA reads)', title = paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                                                                                   ': p-val =',signif(temp$coefficients[2,4],4))) +
  theme_bw()+ theme(legend.position = "none")

p1 | p2 |p3 + theme(legend.position = "right")


dev.off()


pdf("pdf/three/figure6A.pdf", width = 17, height = 5)

par(mfrow=c(1,3))
TF.name='JUN(var.2)'
#boxplot(motif.mat[which(TF==TF.name),]~hiv$group)
temp=summary(lm(sqrt(pATAC)~motif.mat[which(TF==TF.name),]))
temp1 = data.frame( cbind(motif.mat[which(TF==TF.name),],sqrt(pATAC)))
temp2 = paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
              ': p-val =',signif(temp$coefficients[2,4],4))

p1 = ggplot(temp1, aes(motif.mat[which(TF==TF.name),], sqrt(pATAC))) +
  geom_point(aes(color = hiv$group)) +
  geom_smooth(method='lm') +
  scale_color_manual(values = color) +
  labs(x= 'TF motif score', y = 'sqrt(percentage of HIV ATAC reads)', title = paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                                                                                   ': p-val =',signif(temp$coefficients[2,4],4))) +
  theme_bw()+ theme(legend.position = "none")

TF.name='FOS'
#boxplot(motif.mat[which(TF==TF.name),]~hiv$group)
temp=summary(lm(sqrt(pATAC)~motif.mat[which(TF==TF.name),]))
temp1 = data.frame( cbind(motif.mat[which(TF==TF.name),],sqrt(pATAC)))
p2 = ggplot(temp1, aes(motif.mat[which(TF==TF.name),], sqrt(pATAC))) +
  geom_point(aes(color = hiv$group)) +
  geom_smooth(method='lm') +
  scale_color_manual(values = color) +
  labs(x= 'TF motif score', y = 'sqrt(percentage of HIV ATAC reads)', title = paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                                                                                   ': p-val =',signif(temp$coefficients[2,4],4))) +
  theme_bw() + theme(legend.position = "none")


TF.name='FOS::JUN'
temp=summary(lm(sqrt(pATAC)~motif.mat[which(TF==TF.name),]))
temp1 = data.frame( cbind(motif.mat[which(TF==TF.name),],sqrt(pATAC)))
p3 = ggplot(temp1, aes(motif.mat[which(TF==TF.name),], sqrt(pATAC))) +
  geom_point(aes(color = hiv$group)) +
  geom_smooth(method='lm') +
  scale_color_manual(values = color) +
  labs(x= 'TF motif score', y = 'sqrt(percentage of HIV ATAC reads)', title = paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                                                                                   ': p-val =',signif(temp$coefficients[2,4],4))) +
  theme_bw()+ theme(legend.position = "none")

p1 | p2 |p3 + theme(legend.position = "right")


dev.off()


# Get the p-value for each TF
linked.TF=matrix(nrow=length(TF), ncol=4)
linked.TF=as.data.frame(linked.TF)
linked.TF[,1]=TF
linked.TF[,2]=rownames(motif.mat)
colnames(linked.TF)=c('TF','motif','cor','pval')
for(i in 1:length(TF)){
  TF.name=TF[i]
  temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
  linked.TF[i,3]=cor(sqrt(pRNA), motif.mat[which(TF==TF.name),], method='spearman')
  linked.TF[i,4]=signif(temp$coefficients[2,4],4)
}

pdf("pdf/three/figure7.pdf", width = 5, height = 5)
hist(linked.TF[,4], xlab='p-value', main='P-val dist. from linking HIV exp. and TF score')
dev.off()

library(qvalue)
linked.TF=cbind(linked.TF, pval.adj=qvalue(p=linked.TF[,4])$qvalues)

a = -log(linked.TF[,5],10)>-log(0.05,10)
a[a == TRUE] <- "red"
a[a == FALSE] <- "black"
pdf("pdf/three/figure8.pdf", width = 5, height = 5)

plot(linked.TF[,3], -log(linked.TF[,5],10), pch=16, xlab='Spearman correlation', ylab='-log(adj. p-value)', cex=0.6, col=a)

abline(h=-log(0.05,10), lty=2)
temp=linked.TF[order(linked.TF[,5]),]
temp[temp$cor<0,]
dev.off()

write.csv(linked.TF[order(linked.TF[,5]),], file='linkage_analysis_output_wm/linked.TF.csv', row.names = F)



# Get the p-value for each TF
linked.TF=matrix(nrow=length(TF), ncol=4)
linked.TF=as.data.frame(linked.TF)
linked.TF[,1]=TF
linked.TF[,2]=rownames(motif.mat)
colnames(linked.TF)=c('TF','motif','cor','pval')
for(i in 1:length(TF)){
  TF.name=TF[i]
  temp=summary(lm(sqrt(pATAC)~motif.mat[which(TF==TF.name),]))
  linked.TF[i,3]=cor(sqrt(pATAC), motif.mat[which(TF==TF.name),], method='spearman')
  linked.TF[i,4]=signif(temp$coefficients[2,4],4)
}

pdf("pdf/three/figure7.1.pdf", width = 5, height = 5)
hist(linked.TF[,4], xlab='p-value', main='P-val dist. from linking HIV ATAC read counts and TF score')
dev.off()

library(qvalue)
linked.TF=cbind(linked.TF, pval.adj=qvalue(p=linked.TF[,4])$qvalues)

a = -log(linked.TF[,5],10)>-log(0.05,10)
a[a == TRUE] <- "red"
a[a == FALSE] <- "black"
pdf("pdf/three/figure8.1.pdf", width = 5, height = 5)

plot(linked.TF[,3], -log(linked.TF[,5],10), pch=16, xlab='Spearman correlation', ylab='-log(adj. p-value)', col = a,cex=0.6)

abline(h=-log(0.05,10), lty=2)
temp=linked.TF[order(linked.TF[,5]),]
temp[temp$cor<0,]
dev.off()

write.csv(linked.TF[order(linked.TF[,5]),], file='linkage_analysis_output_atac/linked.TF.csv', row.names = F)


linked.TF=matrix(nrow=length(TF), ncol=5)
linked.TF=as.data.frame(linked.TF)
linked.TF[,1]=TF
linked.TF[,2]=rownames(motif.mat)
colnames(linked.TF)=c('TF','motif','RNA.cor','ATAC.cor','pval')
for(i in 1:length(TF)){
  TF.name=TF[i]
  temp=Anova(lm(cbind(sqrt(pATAC),sqrt(pRNA))~motif.mat[which(TF==TF.name),]))
  linked.TF[i,3]=cor(sqrt(pRNA), motif.mat[which(TF==TF.name),], method='spearman')
  linked.TF[i,4]=cor(sqrt(pATAC), motif.mat[which(TF==TF.name),], method='spearman')
  linked.TF[i,5]=as.numeric(sub(".*\\s*(\\d+\\.[0-9e-]+)\\s*[*.]*", "\\1", capture.output(temp)[4]))
}




pdf("pdf/three/figure7.2.pdf", width = 5, height = 5)
hist(linked.TF[,5], xlab='p-value', main='P-val dist. from linking HIV ATAC and RNA read counts jointly to TF score')
dev.off()

library(qvalue)
linked.TF=cbind(linked.TF, pval.adj=qvalue(p=linked.TF[,5])$qvalues)
a = -log(linked.TF[,6],10)>-log(0.05,10)
a[a == TRUE] <- "red"
a[a == FALSE] <- "black"
pdf("pdf/three/figure8.2.pdf", width = 5, height = 5)

plot(linked.TF[,3], -log(linked.TF[,6],10), pch=16, xlab='Spearman correlation', col = a,ylab='-log(adj. p-value)', cex=0.6)
abline(h=-log(0.05,10), lty=2)
temp=linked.TF[order(linked.TF[,6]),]
temp[temp$RNA.cor<0,]

plot(linked.TF[,4], -log(linked.TF[,6],10), pch=16, xlab='Spearman correlation', col = a,ylab='-log(adj. p-value)', cex=0.6)
abline(h=-log(0.05,10), lty=2)

dev.off()

write.csv(linked.TF[order(linked.TF[,5]),], file='linkage_analysis_output_joint/linked.TF.csv', row.names = F)



pdf("pdf/three/figure9.pdf", width = 16, height = 4)
# Plot for top linked TFs
par(mfrow=c(1,4))
for(TF.name in c('FOS::JUND','FOSB::JUNB','FOXP1','SNAI2')){
  temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
  smoothScatter(motif.mat[which(TF==TF.name),],sqrt(pRNA), xlab='TF motif score.', ylab='sqrt(percentage of HIV RNA reads)',
                main=paste(TF.name, '\np-val =',signif(temp$coefficients[2,4],3), 'r =', signif(cor(sqrt(pRNA),motif.mat[which(TF==TF.name),]),3)))
  abline(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]), lty=2, col='black')
  
  # # below is applying regression to only cells with pRNA>0
  # temp=summary(lm(sqrt(pRNA)[pRNA>0]~motif.mat[which(TF==TF.name),][pRNA>0]))
  # smoothScatter(motif.mat[which(TF==TF.name),][pRNA>0],sqrt(pRNA)[pRNA>0], xlab='TF motif score.', ylab='sqrt(percentage of HIV RNA reads)',
  #               main=paste(TF.name, 'p-val =',signif(temp$coefficients[2,4],4)))
  # abline(lm(sqrt(pRNA)[pRNA>0]~motif.mat[which(TF==TF.name),][pRNA>0]), lty=2, col='black')
  
}
par(mfrow=c(1,1))
dev.off()


linked.TF <- read.csv("linkage_analysis_output_atac/linked.TF.csv")
pdf("pdf/three/figure9.1.pdf", width = 12, height = 3)
# Plot for top linked TFs
par(mfrow=c(1,4))
for(TF.name in  c('FOS::JUND','FOSB::JUNB','FOXP1','SNAI2')){
  temp=summary(lm(sqrt(pATAC)~motif.mat[which(TF==TF.name),]))
  smoothScatter(motif.mat[which(TF==TF.name),],sqrt(pATAC), xlab='TF motif score.', ylab='sqrt(percentage of HIV ATAC reads)',
                main=paste(TF.name, '\np-val =',signif(temp$coefficients[2,4],3), 'r =', signif(cor(sqrt(pATAC),motif.mat[which(TF==TF.name),]),3)))
  abline(lm(sqrt(pATAC)~motif.mat[which(TF==TF.name),]), lty=2, col='black')
  
}
par(mfrow=c(1,1))
dev.off()



var.genes=hiv@assays$SCT@var.features

linked.gene=matrix(nrow=length(var.genes), ncol=3)
linked.gene=as.data.frame(linked.gene)
linked.gene[,1]=var.genes
colnames(linked.gene)=c('gene','cor','pval')
for(i in 1:nrow(linked.gene)){
  gene.name=var.genes[i]
  if(i %%200==0) cat(i,'\t')
  temp=summary(lm(sqrt(pRNA)~hiv@assays$SCT@data[gene.name,]))
  linked.gene[i,2]=cor(sqrt(pRNA), hiv@assays$SCT@data[gene.name,], method='spearman')
  linked.gene[i,3]=signif(temp$coefficients[2,4],4)
}


a = -log(linked.gene[,4],10)>-log(0.05,10)
a[a == TRUE] <- "red"
a[a == FALSE] <- "black"

pdf("pdf/three/figure10.pdf", width = 5, height = 5)
hist(linked.gene[,3], xlab='p-value', main='P-val dist. from linking HIV exp. and gene exp.')

library(qvalue)
linked.gene=cbind(linked.gene, pval.adj=qvalue(p=linked.gene[,3])$qvalues)

plot(linked.gene[,2], -log(linked.gene[,4],10), pch=16, xlab='Spearman correlation', ylab='-log(adj. p-value)', cex=0.6,col=a,
     ylim=c(0,20), xlim=c(-0.15, 0.15))
abline(h=-log(0.05,10), lty=2)
dev.off()


linked.gene.output=linked.gene[order(linked.gene[,4]),]
write.table(linked.gene.output, file='linkage_analysis_output_wm/linked.gene.output.txt', row.names = F, col.names = T, sep='\t', quote=FALSE)

write.csv(linked.gene.output, file='linkage_analysis_output_wm/linked.gene.csv', row.names = F)


linked.gene.output




pdf("pdf/three/figure11.pdf",width = 12, height = 3)
# Plot for top linked genes
par(mfrow=c(1,4))
#for(gene.name in c('SORL1','TC2N','MKI67','LINC00861','MALAT1', 'NEAT1')){
for(gene.name in c('MKI67','LTB','CENPF', 'IL7R')){
  temp=summary(lm(sqrt(pRNA)~hiv@assays$SCT@data[gene.name,]))
  smoothScatter(hiv@assays$SCT@data[gene.name,],sqrt(pRNA), xlab='RNA exp.', ylab='sqrt(percentage of HIV RNA reads)',
                main=paste(gene.name, '\np-val =',signif(temp$coefficients[2,4],3), 'r =', signif(cor(sqrt(pRNA),hiv@assays$SCT@data[gene.name,]),3)))
  abline(lm(sqrt(pRNA)~hiv@assays$SCT@data[gene.name,]), lty=2, col='black')
}
par(mfrow=c(1,1))

dev.off()




var.genes=hiv@assays$SCT@var.features

linked.gene=matrix(nrow=length(var.genes), ncol=3)
linked.gene=as.data.frame(linked.gene)
linked.gene[,1]=var.genes
colnames(linked.gene)=c('gene','cor','pval')
for(i in 1:nrow(linked.gene)){
  gene.name=var.genes[i]
  if(i %%200==0) cat(i,'\t')
  temp=summary(lm(sqrt(pATAC)~hiv@assays$SCT@data[gene.name,]))
  linked.gene[i,2]=cor(sqrt(pATAC), hiv@assays$SCT@data[gene.name,], method='spearman')
  linked.gene[i,3]=signif(temp$coefficients[2,4],4)
}

a = -log(linked.gene[,4],10)>-log(0.05,10)
a[a == TRUE] <- "red"
a[a == FALSE] <- "black"
pdf("pdf/three/figure10.1.pdf", width = 5, height = 5)
hist(linked.gene[,3], xlab='p-value', main='P-val dist. from linking HIV ATAC read counts and gene exp.')

library(qvalue)
linked.gene=cbind(linked.gene, pval.adj=qvalue(p=linked.gene[,3])$qvalues)

plot(linked.gene[,2], -log(linked.gene[,4],10), pch=16, xlab='Spearman correlation', ylab='-log(adj. p-value)', cex=0.6,col=a,
     ylim=c(0,20), xlim=c(-0.15, 0.15))
abline(h=-log(0.05,10), lty=2)
dev.off()


linked.gene.output=linked.gene[order(linked.gene[,4]),]
write.table(linked.gene.output, file='linkage_analysis_output_atac/linked.gene.output.txt', row.names = F, col.names = T, sep='\t', quote=FALSE)

write.csv(linked.gene.output, file='linkage_analysis_output_atac/linked.gene.csv', row.names = F)


linked.gene.output

pdf("pdf/three/figure11.1.pdf", width = 12, height = 3)
# Plot for top linked genes
par(mfrow=c(1,4))
#for(gene.name in c('SORL1','TC2N','MKI67','LINC00861','MALAT1', 'NEAT1')){
  for(gene.name in c('MKI67','LTB','CENPF', 'IL7R')){
  temp=summary(lm(sqrt(pATAC)~hiv@assays$SCT@data[gene.name,]))
  smoothScatter(hiv@assays$SCT@data[gene.name,],sqrt(pATAC), xlab='ATAC exp.', ylab='sqrt(percentage of HIV ATAC reads)',
                main=paste(gene.name, '\np-val =',signif(temp$coefficients[2,4],3), 'r =', signif(cor(sqrt(pATAC),hiv@assays$SCT@data[gene.name,]),3)))
  abline(lm(sqrt(pATAC)~hiv@assays$SCT@data[gene.name,]), lty=2, col='black')
}
par(mfrow=c(1,1))

dev.off()

var.genes=hiv@assays$SCT@var.features
linked.gene=matrix(nrow=length(var.genes), ncol=4)
linked.gene=as.data.frame(linked.gene)
linked.gene[,1]=var.genes
colnames(linked.gene)=c('gene','RNA.cor','ATAC.cor','pval')
for(i in 1:nrow(linked.gene)){
  gene.name=var.genes[i]
  if(i %%200==0) cat(i,'\t')
  temp=Anova(lm(cbind(sqrt(pATAC),sqrt(pRNA))~hiv@assays$SCT@data[gene.name,]))
  linked.gene[i,3]=cor(sqrt(pATAC), hiv@assays$SCT@data[gene.name,], method='spearman')
  linked.gene[i,2]=cor(sqrt(pRNA), hiv@assays$SCT@data[gene.name,], method='spearman')
  linked.gene[i,4]=as.numeric(sub(".*\\s*(\\d+\\.[0-9e-]+)\\s*[*.]*", "\\1", capture.output(temp)[4]))
}

pdf("pdf/three/figure10.2.pdf", width = 5, height = 5)
hist(linked.gene[,4], xlab='p-value', main='P-val dist. from linking HIV ATAC and RNA read counts jointly and gene exp.')

a = -log(linked.gene[,5],10)>-log(0.05,10)
a[a == TRUE] <- "red"
a[a == FALSE] <- "black"
library(qvalue)
linked.gene=cbind(linked.gene, pval.adj=qvalue(p=linked.gene[,4])$qvalues)

plot(linked.gene[,2], -log(linked.gene[,5],10), pch=16, xlab='Spearman correlation', ylab='-log(adj. p-value)', cex=0.6,
     ylim=c(0,20), xlim=c(-0.15, 0.15), col=a)
abline(h=-log(0.05,10), lty=2)

plot(linked.gene[,3], -log(linked.gene[,5],10), pch=16, xlab='Spearman correlation', ylab='-log(adj. p-value)', cex=0.6,
     ylim=c(0,20), xlim=c(-0.15, 0.15), col=a)
abline(h=-log(0.05,10), lty=2)
dev.off()


linked.gene.output=linked.gene[order(linked.gene[,5]),]
write.table(linked.gene.output, file='linkage_analysis_output_joint/linked.gene.output.txt', row.names = F, col.names = T, sep='\t', quote=FALSE)

write.csv(linked.gene.output, file='linkage_analysis_output_joint/linked.gene.csv', row.names = F)


linked.gene.output


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
  linked.TF.gene[i,2]=cor(sqrt(pRNA), hiv@assays$SCT@data[gene.name,], method='spearman')
  linked.TF.gene[i,3]=signif(temp$coefficients[2,4],4)
}
pdf("pdf/three/figure12.pdf", width = 10, height = 8)
hist(linked.TF.gene[,3], xlab='p-value', main='P-val dist. from linking HIV exp. and TF gene exp.')

library(qvalue)
linked.TF.gene=cbind(linked.TF.gene, pval.adj=qvalue(p=linked.TF.gene[,3])$qvalues)

plot(linked.TF.gene[,2], -log(linked.TF.gene[,4],10), pch=16, xlab='Spearman correlation', ylab='-log(adj. p-value)', cex=0.6,
     ylim=c(0,4), xlim=c(-0.15, 0.15))
abline(h=-log(0.05,10), lty=2)

dev.off()

linked.TF.gene.output=linked.TF.gene[order(linked.TF.gene[,4]),]
colnames(linked.TF.gene.output)[1]='TF.gene'
write.table(linked.TF.gene.output, file='linkage_analysis_output_wm/linked.TF.gene.output.txt', row.names = F, col.names = T, sep='\t', quote=FALSE)

write.csv(linked.TF.gene, file='linkage_analysis_output_wm/linked.TF.gene.csv', row.names = F)
linked.TF.gene.output[which(linked.TF.gene.output$cor <0),]
# Plot for top linked TF genes

pdf("pdf/three/figure13.pdf", width = 10, height = 8)
par(mfrow=c(2,2))
for(gene.name in c('TCF7','LEF1','NR1D2','KLF2')){
  temp=summary(lm(sqrt(pRNA)~hiv@assays$SCT@data[gene.name,]))
  smoothScatter(hiv@assays$SCT@data[gene.name,],sqrt(pRNA), xlab='TF RNA exp.', ylab='sqrt(percentage of HIV RNA reads)',
                main=paste(gene.name, 'p-val =',signif(temp$coefficients[2,4],4)))
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
  temp=summary(lm(sqrt(pATAC)~hiv@assays$SCT@data[gene.name,]))
  linked.TF.gene[i,2]=cor(sqrt(pATAC), hiv@assays$SCT@data[gene.name,], method='spearman')
  linked.TF.gene[i,3]=signif(temp$coefficients[2,4],4)
}


pdf("pdf/three/figure12.1.pdf", width = 5, height = 5)
hist(linked.TF.gene[,3], xlab='p-value', main='P-val dist. from linking HIV ATAC read counts and TF gene exp.')

library(qvalue)
linked.TF.gene=cbind(linked.TF.gene, pval.adj=qvalue(p=linked.TF.gene[,3])$qvalues)

plot(linked.TF.gene[,2], -log(linked.TF.gene[,4],10), pch=16, xlab='Spearman correlation', ylab='-log(adj. p-value)', cex=0.6,
     ylim=c(0,4), xlim=c(-0.15, 0.15))
abline(h=-log(0.05,10), lty=2)

dev.off()

linked.TF.gene.output=linked.TF.gene[order(linked.TF.gene[,4]),]
colnames(linked.TF.gene.output)[1]='TF.gene'
write.table(linked.TF.gene.output, file='linkage_analysis_output_atac/linked.TF.gene.output.txt', row.names = F, col.names = T, sep='\t', quote=FALSE)

write.csv(linked.TF.gene, file='linkage_analysis_output_atac/linked.TF.gene.csv', row.names = F)
linked.TF.gene.output[which(linked.TF.gene.output$cor <0),]
# Plot for top linked TF genes

pdf("pdf/three/figure13.1.pdf", width = 12, height = 3)
par(mfrow=c(1,4))
for(gene.name in c('TCF7','LEF1','NR1D2','KLF2')){
  temp=summary(lm(sqrt(pATAC)~hiv@assays$SCT@data[gene.name,]))
  smoothScatter(hiv@assays$SCT@data[gene.name,],sqrt(pATAC), xlab='TF ATAC exp.', ylab='sqrt(percentage of HIV ATAC reads)',
                main=paste(gene.name, 'p-val =',signif(temp$coefficients[2,4],4)))
  abline(lm(sqrt(pATAC)~hiv@assays$SCT@data[gene.name,]), lty=2, col='black')
}
par(mfrow=c(1,1))
dev.off()




# Get the p-value for each TF gene expression
TF.genes=TF[TF %in% rownames(hiv@assays$SCT)]

linked.TF.gene=matrix(nrow=length(TF.genes), ncol=4)
linked.TF.gene=as.data.frame(linked.TF.gene)
linked.TF.gene[,1]=TF.genes
colnames(linked.TF.gene)=c('gene','RNA.cor','ATAC.cor','pval')
for(i in 1:nrow(linked.TF.gene)){
  gene.name=TF.genes[i]
  if(i %%200==0) cat(i,'\t')
  temp=Anova(lm(cbind(sqrt(pATAC),sqrt(pRNA))~hiv@assays$SCT@data[gene.name,]))
  linked.TF.gene[i,2]=cor(sqrt(pRNA), hiv@assays$SCT@data[gene.name,], method='spearman')
  linked.TF.gene[i,3]=cor(sqrt(pATAC), hiv@assays$SCT@data[gene.name,], method='spearman')
  linked.TF.gene[i,4]=as.numeric(sub(".*\\s*(\\d+\\.[0-9e-]+)\\s*[*.]*", "\\1", capture.output(temp)[4]))
}


pdf("pdf/three/figure12.2.pdf", width = 10, height = 8)
hist(linked.TF.gene[,4], xlab='p-value', main='P-val dist. from linking HIV ATAC and RNA read counts jointly to TF gene exp.')

library(qvalue)
linked.TF.gene=cbind(linked.TF.gene, pval.adj=qvalue(p=linked.TF.gene[,4])$qvalues)

plot(linked.TF.gene[,2], -log(linked.TF.gene[,5],10), pch=16, xlab='Spearman correlation', ylab='-log(adj. p-value)', cex=0.6,
     ylim=c(0,4), xlim=c(-0.15, 0.15))
abline(h=-log(0.05,10), lty=2)

dev.off()

linked.TF.gene.output=linked.TF.gene[order(linked.TF.gene[,5]),]
colnames(linked.TF.gene.output)[1]='TF.gene'
write.table(linked.TF.gene.output, file='linkage_analysis_output_joint/linked.TF.gene.output.txt', row.names = F, col.names = T, sep='\t', quote=FALSE)

write.csv(linked.TF.gene, file='linkage_analysis_output_joint/linked.TF.gene.csv', row.names = F)
linked.TF.gene.output[which(linked.TF.gene.output$RNA.cor <0),]



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
  linked.peak[i,2]=cor(sqrt(pRNA), hiv@assays$ATAC@data[peak.index,], method='spearman')
  linked.peak[i,3]=signif(temp$coefficients[2,4],4)
}

pdf("pdf/three/figure14.pdf", width = 5, height = 5)
a = -log(linked.peak[,3],10)>-log(0.01,10)
a[a == TRUE] <- "red"
a[a == FALSE] <- "black"
hist(linked.peak[,3], xlab='p-value', main='P-val dist. from linking HIV exp. and peak access.',breaks = 10)



#  We will just use nominal p-values
plot(linked.peak[,2], -log(linked.peak[,3],10), pch=16, xlab='Spearman correlation', ylab='-log(p-value)', cex=0.6, col=a,
     ylim=c(0,4), xlim=c(-0.1, 0.1))
abline(h=-log(0.01,10), lty=2)
dev.off()

linked.peak.output=linked.peak[order(linked.peak[,3]),] # only output significant ones
write.table(linked.peak.output, file='linkage_analysis_output_wm/linked.peak.output.txt', row.names = F, col.names = T, sep='\t', quote=FALSE)

write.csv(linked.peak, file='linkage_analysis_output_wm/linked.peak.csv', row.names = F)


# Get the p-value for each peak
var.peaks=hiv@assays$ATAC@var.features[1:20000]

linked.peak=matrix(nrow=length(var.peaks), ncol=3)
linked.peak=as.data.frame(linked.peak)
linked.peak[,1]=var.peaks
colnames(linked.peak)=c('peak','cor','pval')
for(i in 1:nrow(linked.peak)){
  peak.index=which(rownames(hiv@assays$ATAC)==var.peaks[i])
  if(i %%400==0) cat(i,'\t')
  temp=summary(lm(sqrt(pATAC)~hiv@assays$ATAC@data[peak.index,]))
  linked.peak[i,2]=cor(sqrt(pATAC), hiv@assays$ATAC@data[peak.index,], method='spearman')
  linked.peak[i,3]=signif(temp$coefficients[2,4],4)
}

pdf("pdf/three/figure14.1.pdf", width = 5, height = 5)
a = -log(linked.peak[,3],10)>-log(0.01,10)
a[a == TRUE] <- "red"
a[a == FALSE] <- "black"
hist(linked.peak[,3], xlab='p-value', main='P-val dist. from linking HIV ATAC read counts and peak access.',breaks =  10)


#  We will just use nominal p-values
plot(linked.peak[,2], -log(linked.peak[,3],10), pch=16, xlab='Spearman correlation', ylab='-log(p-value)', cex=0.6,col=a,
     ylim=c(0,4), xlim=c(-0.1, 0.1))
abline(h=-log(0.01,10), lty=2)
dev.off()

linked.peak.output=linked.peak[order(linked.peak[,3]),] # only output significant ones
write.table(linked.peak.output, file='linkage_analysis_output_atac/linked.peak.output.txt', row.names = F, col.names = T, sep='\t', quote=FALSE)

write.csv(linked.peak, file='linkage_analysis_output_atac/linked.peak.csv', row.names = F)

# Get the p-value for each peak
var.peaks=hiv@assays$ATAC@var.features[1:20000]

linked.peak=matrix(nrow=length(var.peaks), ncol=4)
linked.peak=as.data.frame(linked.peak)
linked.peak[,1]=var.peaks
colnames(linked.peak)=c('peak','RNA.cor','ATAC.cor','pval')
for(i in 1:nrow(linked.peak)){
  peak.index=which(rownames(hiv@assays$ATAC)==var.peaks[i])
  if(i %%400==0) cat(i,'\t')
  temp=Anova(lm(cbind(sqrt(pATAC),sqrt(pRNA))~hiv@assays$ATAC@data[peak.index,]))
  linked.peak[i,2]=cor(sqrt(pRNA), hiv@assays$ATAC@data[peak.index,], method='spearman')
  linked.peak[i,3]=cor(sqrt(pATAC), hiv@assays$ATAC@data[peak.index,], method='spearman')
  linked.peak[i,4]=as.numeric(sub(".*\\s*(\\d+\\.[0-9e-]+)\\s*[*.]*", "\\1", capture.output(temp)[4]))
}

pdf("pdf/three/figure14.2.pdf", width = 5, height = 5)
hist(linked.peak[,4], xlab='p-value', main='P-val dist. from linking HIV ATAC and RNA read counts jointly and peak access.',breaks = 10)

a = -log(linked.peak[,4],10)>-log(0.01,10)
a[a == TRUE] <- "red"
a[a == FALSE] <- "black"

#  We will just use nominal p-values
plot(linked.peak[,2], -log(linked.peak[,4],10), pch=16, xlab='Spearman correlation', ylab='-log(p-value)', cex=0.6,
     ylim=c(0,4), xlim=c(-0.1, 0.1), col=a)
abline(h=-log(0.01,10), lty=2)
plot(linked.peak[,3], -log(linked.peak[,4],10), pch=16, xlab='Spearman correlation', ylab='-log(p-value)', cex=0.6,
     ylim=c(0,4), xlim=c(-0.1, 0.1), col=a)
abline(h=-log(0.01,10), lty=2)
dev.off()

linked.peak.output=linked.peak[order(linked.peak[,4]),] # only output significant ones
write.table(linked.peak.output, file='linkage_analysis_output_joint/linked.peak.output.txt', row.names = F, col.names = T, sep='\t', quote=FALSE)

write.csv(linked.peak, file='linkage_analysis_output_joint/linked.peak.csv', row.names = F)


linked.peak.output = read.table("linkage_analysis_output_wm/linked.peak.output.txt", header = T)
#I'm no sure what this line is for, it keeps producing null 
load("hiv_chromvar.rda")
seqinfo(hiv)=seqinfo(hiv)[seqnames(hiv)[1:25]]
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38)='UCSC'
# add motif information
hiv <- AddMotifs(
  object = hiv,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)


linked.peak.output = read.table("linkage_analysis_output_wm/linked.peak.output.txt", header = T)
sig.link.peaks= linked.peak.output[linked.peak.output[,3]<=0.05 ,1]
sig.link.peaks=sig.link.peaks[!grepl('HIV', sig.link.peaks)] 
sig.link.peaks=sig.link.peaks[!grepl('MT', sig.link.peaks)]


enriched.motifs <- FindMotifs(
  object = hiv,
  features =sig.link.peaks
)


enriched.motifs
write.csv(enriched.motifs, file='linkage_analysis_output_wm/enriched.motif.csv', row.names = F)


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
write.csv(pos.enriched.motifs, file='linkage_analysis_output_wm/pos.enriched.motif.csv', row.names = F)
pdf("pdf/three/figure16.pdf", width = 15, height = 5)
p=MotifPlot(
  object = hiv,
  motifs = head(rownames(pos.enriched.motifs)),
  ncol = 6
)
p
dev.off()

#Negatively linked peaks
sig.neg.link.peaks= linked.peak.output[linked.peak.output$pval<=0.01 & linked.peak.output$cor<0,1]
sig.neg.link.peaks=sig.neg.link.peaks[!grepl('HIV', sig.neg.link.peaks)]
sig.neg.link.peaks=sig.neg.link.peaks[!grepl('MT', sig.neg.link.peaks)]

set.seed(1)
neg.enriched.motifs <- FindMotifs(
  object = hiv,
  features =sig.neg.link.peaks
)
neg.enriched.motifs
write.csv(neg.enriched.motifs, file='linkage_analysis_output_wm/neg.enriched.motif.csv', row.names = F)

pdf("pdf/three/figure17.pdf", width = 15, height = 5)
p=MotifPlot(
  object = hiv,
  motifs = head(rownames(neg.enriched.motifs)),
  ncol = 6
)
p
dev.off()

pdf("pdf/three/figure15.pdf", width = 8, height = 12)
MotifPlot(
  object = hiv,
  motifs = head(rownames(enriched.motifs)),
  ncol = 2
)

dev.off()


linked.peak.output = read.table("linkage_analysis_output_atac/linked.peak.output.txt", header = T)
sig.link.peaks= linked.peak.output[linked.peak.output[,3]<=0.05 ,1]
sig.link.peaks=sig.link.peaks[!grepl('HIV', sig.link.peaks)] 
sig.link.peaks=sig.link.peaks[!grepl('MT', sig.link.peaks)]


enriched.motifs <- FindMotifs(
  object = hiv,
  features =sig.link.peaks
)


enriched.motifs
write.csv(enriched.motifs, file='linkage_analysis_output_atac/enriched.motif.csv', row.names = F)

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
write.csv(pos.enriched.motifs, file='linkage_analysis_output_atac/pos.enriched.motif.csv', row.names = F)
pdf("pdf/three/figure16.1.pdf", width = 15, height = 5)
p=MotifPlot(
  object = hiv,
  motifs = head(rownames(pos.enriched.motifs)),
  ncol = 6
)
p
dev.off()

#Negatively linked peaks
sig.neg.link.peaks= linked.peak.output[linked.peak.output$pval<=0.01 & linked.peak.output$cor<0,1]
sig.neg.link.peaks=sig.neg.link.peaks[!grepl('HIV', sig.neg.link.peaks)]
sig.neg.link.peaks=sig.neg.link.peaks[!grepl('MT', sig.neg.link.peaks)]

set.seed(1)
neg.enriched.motifs <- FindMotifs(
  object = hiv,
  features =sig.neg.link.peaks
)
neg.enriched.motifs
write.csv(neg.enriched.motifs, file='linkage_analysis_output_atac/neg.enriched.motif.csv', row.names = F)

pdf("pdf/three/figure17.1.pdf", width = 15, height = 5)
p=MotifPlot(
  object = hiv,
  motifs = head(rownames(neg.enriched.motifs)),
  ncol = 6
)
p
dev.off()

pdf("pdf/three/figure15.1.pdf", width = 8, height = 12)
MotifPlot(
  object = hiv,
  motifs = head(rownames(enriched.motifs)),
  ncol = 2
)

dev.off()


linked.peak.output = read.table("linkage_analysis_output_joint/linked.peak.output.txt", header = T)

sig.link.peaks= linked.peak.output[linked.peak.output[,4]<=0.05 ,1]
sig.link.peaks=sig.link.peaks[!grepl('HIV', sig.link.peaks)] 
sig.link.peaks=sig.link.peaks[!grepl('MT', sig.link.peaks)]


enriched.motifs <- FindMotifs(
  object = hiv,
  features =sig.link.peaks
)


enriched.motifs
write.csv(enriched.motifs, file='linkage_analysis_output_joint/enriched.motif.csv', row.names = F)

# Positively linked peaks
sig.pos.link.peaks= linked.peak.output[linked.peak.output$pval<=0.01 & linked.peak.output$RNA.cor>0 & linked.peak.output$ATAC.cor>0,1]
sig.pos.link.peaks=sig.pos.link.peaks[!grepl('HIV', sig.pos.link.peaks)]
sig.pos.link.peaks=sig.pos.link.peaks[!grepl('MT', sig.pos.link.peaks)]

set.seed(1)
pos.enriched.motifs <- FindMotifs(
  object = hiv,
  features =sig.pos.link.peaks
)
pos.enriched.motifs
write.csv(pos.enriched.motifs, file='linkage_analysis_output_joint/pos.enriched.motif.csv', row.names = F)
pdf("pdf/three/figure16.2.pdf", width = 15, height = 5)
p=MotifPlot(
  object = hiv,
  motifs = head(rownames(pos.enriched.motifs)),
  ncol = 6
)
p
dev.off()

#Negatively linked peaks
sig.neg.link.peaks= linked.peak.output[linked.peak.output$pval<=0.01 & linked.peak.output$RNA.cor<0  & linked.peak.output$ATAC.cor<0,1]
sig.neg.link.peaks=sig.neg.link.peaks[!grepl('HIV', sig.neg.link.peaks)]
sig.neg.link.peaks=sig.neg.link.peaks[!grepl('MT', sig.neg.link.peaks)]

set.seed(1)
neg.enriched.motifs <- FindMotifs(
  object = hiv,
  features =sig.neg.link.peaks
)
neg.enriched.motifs
write.csv(neg.enriched.motifs, file='linkage_analysis_output_joint/neg.enriched.motif.csv', row.names = F)

pdf("pdf/three/figure17.2.pdf", width = 15, height = 5)
p=MotifPlot(
  object = hiv,
  motifs = head(rownames(neg.enriched.motifs)),
  ncol = 6
)
p
dev.off()

pdf("pdf/three/figure15.2.pdf", width = 8, height = 12)
MotifPlot(
  object = hiv,
  motifs = head(rownames(enriched.motifs)),
  ncol = 2
)

dev.off()



linked.TF.rna <- read.csv('linkage_analysis_output_wm/linked.TF.csv')
linked.TF.atac<- read.csv('linkage_analysis_output_atac/linked.TF.csv')
linked.TF.joint <- read.csv('linkage_analysis_output_joint/linked.TF.csv')

library(ggvenn)
x=list(RNA=linked.TF.rna[linked.TF.rna$pval.adj <= 0.1,2],
       ATAC = linked.TF.atac[linked.TF.atac$pval.adj <= 0.1,2],
       Joint = linked.TF.joint[linked.TF.joint$pval.adj <= 0.1,2])

p1=ggvenn(x, 
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)
p1

linked.gene.rna <- read.csv('linkage_analysis_output_wm/linked.gene.csv')
linked.gene.atac<- read.csv('linkage_analysis_output_atac/linked.gene.csv')
linked.gene.joint <- read.csv('linkage_analysis_output_joint/linked.gene.csv')
library(ggvenn)
x=list(RNA=linked.gene.rna[linked.gene.rna$pval.adj <= 0.1,1],
       ATAC = linked.gene.atac[linked.gene.atac$pval.adj <= 0.1,1],
       Joint = linked.gene.joint[linked.gene.joint$pval.adj <= 0.1,1])

p2=ggvenn(x, 
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)
p2

enriched.rna <- read.csv('linkage_analysis_output_wm/enriched.motif.csv')
enriched.atac <- read.csv('linkage_analysis_output_atac/enriched.motif.csv')
enriched.joint <-  read.csv('linkage_analysis_output_joint/enriched.motif.csv')

library(ggvenn)
x=list(RNA=enriched.rna[enriched.rna$p.adjust <= 0.1,1],
       ATAC = enriched.atac[enriched.atac$p.adjust <= 0.1,1],
       Joint = enriched.joint[enriched.joint$p.adjust <= 0.1,1])

p4=ggvenn(x, 
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)
p4
pdf("pdf/three/venn.enriched.pdf", width = 10, height = 10)
p4
dev.off()

linked.peak.rna <- read.table('linkage_analysis_output_wm/linked.peak.output.txt', header = T)
linked.peak.atac<- read.table('linkage_analysis_output_atac/linked.peak.output.txt', header = T)
linked.peak.joint <- read.table('linkage_analysis_output_joint/linked.peak.output.txt', header = T)
library(ggvenn)
x=list(RNA=linked.peak.rna[linked.peak.rna$pval <= 0.01,1],
       ATAC = linked.peak.atac[linked.peak.atac$pval <= 0.01,1],
       Joint = linked.peak.joint[linked.peak.joint$pval <= 0.01,1])

p3=ggvenn(x, 
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)
p3

pdf("pdf/three/venn.pdf", width = 10, height = 10)
p1
p2
p3
p4
dev.off()






library(ggvenn)

linked.TF.rna  <- read.csv('linkage_analysis_output_wm//linked.TF.csv')
linked.TF.rna2 <- read.csv('../../HIV2/linkage_analysis_output_wm//linked.TF.csv')

x=list(Donor1_Joint = linked.TF.rna[linked.TF.rna$pval <= 0.01,1],
       Donor2_Joint = linked.TF.rna2[linked.TF.rna2$pval <= 0.01,1])

p1=ggvenn(x,
          fill_color = c("#0073C2FF", "#EFC000FF"),
          stroke_size = 0.5, set_name_size = 2, show_percentage = FALSE
)
p1


linked.gene.rna  <- read.csv('linkage_analysis_output_wm//linked.gene.csv')
linked.gene.rna2 <- read.csv('../../HIV2/linkage_analysis_output_wm/linked.gene.csv')

x=list(Donor1_Joint = linked.gene.rna[linked.gene.rna$pval <= 0.01,1],
       Donor2_Joint = linked.gene.rna2[linked.gene.rna2$pval <= 0.01,1])

p2=ggvenn(x,
          fill_color = c("#0073C2FF", "#EFC000FF"),
          stroke_size = 0.5, set_name_size = 2, show_percentage = FALSE
)
p2

linked.peak.rna  <- read.csv('linkage_analysis_output_wm/enriched.motif.csv')
linked.peak.rna2 <- read.csv('../../HIV2/linkage_analysis_output_wm/enriched.motif.csv')


x=list(Donor1_Joint = linked.peak.rna[linked.peak.rna$pvalue <= 0.01,1],
       Donor2_Joint = linked.peak.rna2[linked.peak.rna2$pvalue <= 0.01,1])

p3=ggvenn(x,
          fill_color = c("#0073C2FF", "#EFC000FF"),
          stroke_size = 0.5, set_name_size = 2, show_percentage = FALSE
)
p3





pdf("pdf/three/venn2.pdf", width = 13, height = 10)
p1
p2
p3

dev.off()



gene1 = "LEF1"
linked.gene.atac[linked.gene.atac$gene==gene1,]
linked.gene.rna[linked.gene.rna$gene==gene1,]
linked.gene.joint[linked.gene.joint$gene==gene1,]

linked.gene.exp.atac[linked.gene.exp.atac$gene==gene1,]
linked.gene.exp.rna[linked.gene.exp.rna$gene==gene1,]
linked.gene.exp.joint[linked.gene.exp.joint$gene==gene1,]


gene2 = "ZNF135"
linked.TF.atac[linked.TF.atac$TF==gene1|linked.TF.atac$TF==gene2,]
linked.TF.rna[linked.TF.rna$TF==gene1|linked.TF.rna$TF==gene2,]
linked.TF.joint[linked.TF.joint$TF==gene1|linked.TF.joint$TF==gene2,]




linked.gene.atac[substring(linked.gene.atac$gene,1,3) =="ZNF",]





pdf("pdf/three/featureplot.pdf", width = 10, height = 10)
FeaturePlot(hiv, features = "pRNA.sqrt")
DefaultAssay(hiv) <- "ATAC"
dev.off()


pdf("pdf/three/gata3.pdf", width = 10, height = 10)
par(mfrow=c(2,2))
TF.name='GATA3'
#boxplot(motif.mat[which(TF==TF.name),]~hiv$group)
temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
smoothScatter(motif.mat[which(TF==TF.name),], sqrt(pRNA),col = color[colnum],  xlab='TF motif score',
              pch = 20, cex = 1, ylab=paste('sqrt(percentage of HIV RNA reads)'),
              main=paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                         ': p-val =',signif(temp$coefficients[2,4],4)))
abline(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]), col='black',lty=2)

#boxplot(motif.mat[which(TF==TF.name),]~hiv$group)
temp=summary(lm(sqrt(pATAC)~motif.mat[which(TF==TF.name),]))
smoothScatter(motif.mat[which(TF==TF.name),], sqrt(pATAC),col = color[colnum],  xlab='TF motif score',
              pch = 20, cex = 1, ylab=paste('sqrt(percentage of HIV ATAC reads)'),
              main=paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                         ': p-val =',signif(temp$coefficients[2,4],4)))
abline(lm(sqrt(pATAC)~motif.mat[which(TF==TF.name),]), col='black',lty=2)


TF.name='FOXP1'
#boxplot(motif.mat[which(TF==TF.name),]~hiv$group)
temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
smoothScatter(motif.mat[which(TF==TF.name),], sqrt(pRNA),col = color[colnum],  xlab='TF motif score',
              pch = 20, cex = 1,  ylab=paste('sqrt(percentage of HIV RNA reads)'),
              main=paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                         ': p-val =',signif(temp$coefficients[2,4],4)))
abline(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]), col='black',lty=2)

temp=summary(lm(sqrt(pATAC)~motif.mat[which(TF==TF.name),]))
smoothScatter(motif.mat[which(TF==TF.name),], sqrt(pATAC),col = color[colnum],  xlab='TF motif score',
              pch = 20, cex = 1, ylab=paste('sqrt(percentage of HIV ATAC reads)'),
              main=paste('TF',TF.name,', motif', rownames(motif.mat)[which(TF==TF.name)],
                         ': p-val =',signif(temp$coefficients[2,4],4)))
abline(lm(sqrt(pATAC)~motif.mat[which(TF==TF.name),]), col='black',lty=2)

dev.off()








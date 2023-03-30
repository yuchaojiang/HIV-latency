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

load("/Users/mwen/Documents/Dissertation/HIV/HIV2/processed.hiv.rda")
#load("/Users/mwen/Documents/Dissertation/HIV/HIV1/downsample/hiv.rda")
DefaultAssay(hiv)='SCT'
table(hiv$group)


pdf("pdf/four/figurebee.pdf", width = 5, height = 5)
par(mfrow=c(1,3))
  i = "HIV"
  temp = as.data.frame(SCT[rownames(SCT)==i,])
  temp$Group = hiv$group       
  names(temp) <-c("Normalized Counts", "Group")
  p<-ggplot(temp, aes(x=Group, y=`Normalized Counts`, color=Group)) +
    geom_violin(trim=FALSE) +labs(title=i) +NoLegend()
  p
  i = "LEF1"
  temp = as.data.frame(SCT[rownames(SCT)==i,])
  temp$Group = hiv$group       
  names(temp) <-c("Normalized Counts", "Group")
  p<-ggplot(temp, aes(x=Group, y=`Normalized Counts`, color=Group)) +
    geom_violin(trim=FALSE)+labs(title=i)+NoLegend()
  p
  i = "FRYL"
  temp = as.data.frame(SCT[rownames(SCT)==i,])
  temp$Group = hiv$group       
  names(temp) <-c("Normalized Counts", "Group")
  p<-ggplot(temp, aes(x=Group, y=`Normalized Counts`, color=Group)) +
    geom_violin(trim=FALSE)+labs(title=i) +NoLegend()
  p

dev.off()


pdf("pdf/four/figurebee2.pdf", width = 5, height = 5)
par(mfrow=c(1,3))
i = "RPL41"
temp = as.data.frame(SCT[rownames(SCT)==i,])
temp$Group = hiv$group       
names(temp) <-c("Normalized Counts", "Group")
p<-ggplot(temp, aes(x=Group, y=`Normalized Counts`, color=Group)) +
  geom_violin(trim=FALSE) +labs(title=i) +NoLegend()
p
i = "RPL37A"
temp = as.data.frame(SCT[rownames(SCT)==i,])
temp$Group = hiv$group       
names(temp) <-c("Normalized Counts", "Group")
p<-ggplot(temp, aes(x=Group, y=`Normalized Counts`, color=Group)) +
  geom_violin(trim=FALSE)+labs(title=i)+NoLegend()
p
i = "RPL30"
temp = as.data.frame(SCT[rownames(SCT)==i,])
temp$Group = hiv$group       
names(temp) <-c("Normalized Counts", "Group")
p<-ggplot(temp, aes(x=Group, y=`Normalized Counts`, color=Group)) +
  geom_violin(trim=FALSE)+labs(title=i) +NoLegend()
p

dev.off()

SAHA.RNA.markers <- FindMarkers(hiv, ident.1 = 'SAHA', ident.2='DMSO',min.pct = 0.1, logfc.threshold = 0.1)
SAHA.RNA.markers=SAHA.RNA.markers[SAHA.RNA.markers$p_val_adj<=0.05,]
SAHA.RNA.markers.down = SAHA.RNA.markers[SAHA.RNA.markers$avg_log2FC<0,]
SAHA.RNA.markers.up = SAHA.RNA.markers[SAHA.RNA.markers$avg_log2FC>0,]

write.csv(SAHA.RNA.markers.down, file='SAHA.RNA.markers.down.csv', row.names = T, quote = F)
write.csv(SAHA.RNA.markers.up, file='AHA.RNA.markers.up.csv', row.names = T, quote = F)

Pros.RNA.markers <- FindMarkers(hiv, ident.1 = 'Prostratin', ident.2='DMSO',min.pct = 0.1, logfc.threshold = 0.1)
Pros.RNA.markers=Pros.RNA.markers[Pros.RNA.markers$p_val_adj<=0.05,]
Pros.RNA.markers.down = Pros.RNA.markers[Pros.RNA.markers$avg_log2FC<0,]
Pros.RNA.markers.up = Pros.RNA.markers[Pros.RNA.markers$avg_log2FC>0,]

write.csv(Pros.RNA.markers.down, file='Pros.RNA.markers.down.csv', row.names = T, quote = F)
write.csv(Pros.RNA.markers.up, file='Pros.RNA.markers.up.csv', row.names = T, quote = F)

iBET.RNA.markers <- FindMarkers(hiv, ident.1 = 'iBET151', ident.2='DMSO',min.pct = 0.1, logfc.threshold = 0.1)
iBET.RNA.markers=iBET.RNA.markers[iBET.RNA.markers$p_val_adj<=0.05,]
iBET.RNA.markers.down = iBET.RNA.markers[iBET.RNA.markers$avg_log2FC<0,]
iBET.RNA.markers.up = iBET.RNA.markers[iBET.RNA.markers$avg_log2FC>0,]

write.csv(iBET.RNA.markers.down, file='iBET.RNA.markers.down.csv', row.names = T, quote = F)
write.csv(iBET.RNA.markers.up, file='iBET.RNA.markers.up.csv', row.names = T, quote = F)



library(ggvenn)
x=list(SAHA.down=rownames(SAHA.RNA.markers.down),
       Prostratin.down=rownames(Pros.RNA.markers.down),
       iBET151.down=rownames(iBET.RNA.markers.down))
p1=ggvenn(x, 
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)

x=list(SAHA.up=rownames(SAHA.RNA.markers.up),
       Prostratin.up=rownames(Pros.RNA.markers.up),
       iBET151.up=rownames(iBET.RNA.markers.up))
p2=ggvenn(x, 
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)

p=p1+p2+plot_layout(ncol=1)
p
ggsave(p, file="pdf/four/figure1.pdf", width=3, height=6)



up.genes=unique(c(rownames(SAHA.RNA.markers.up)[1:15],
                  rownames(Pros.RNA.markers.up)[1:15],
                  rownames(iBET.RNA.markers.up)[1:15]))
p.up=DoHeatmap(hiv, features = up.genes, group.by = 'group')
ggsave("pdf/four/figure1.1.pdf", plot=p.up, width = 5, height=8)

down.genes=unique(c(rownames(SAHA.RNA.markers.down)[1:15],
                    rownames(Pros.RNA.markers.down)[1:15],
                    rownames(iBET.RNA.markers.down)[1:15]))
p.down=DoHeatmap(hiv, features = down.genes, group.by = 'group')
ggsave("pdf/four/figure1.2.pdf", plot=p.down, width = 5, height=8)



pdf("pdf/four/figure1.2.pdf", width = 7, height = 15)
a= rownames(SAHA.RNA.markers.up)
b= rownames(Pros.RNA.markers.up)
c= rownames(iBET.RNA.markers.up)
#up.genes = unique(c(intersect(a,b),intersect(b,c),intersect(c,a)))
up.genes = intersect(intersect(a,b),c)[1:50]


write.table(a,"SAHA.RNA.markers.up.txt",row.names=FALSE,sep="\t", quote = FALSE)
write.table(b,"Pros.RNA.markers.up.txt",row.names=FALSE,sep="\t", quote = FALSE)
write.table(c,"iBET.RNA.markers.up.txt",row.names=FALSE,sep="\t", quote = FALSE)


a= rownames(SAHA.RNA.markers.down)
b= rownames(Pros.RNA.markers.down)
c= rownames(iBET.RNA.markers.down)
#down.genes = unique(c(intersect(a,b),intersect(b,c),intersect(c,a)))
down.genes = intersect(intersect(a,b),c) [1:50]

write.table(a,"SAHA.RNA.markers.down.txt",row.names=FALSE,sep="\t", quote = FALSE)
write.table(b,"Pros.RNA.markers.down.txt",row.names=FALSE,sep="\t", quote = FALSE)
write.table(c,"iBET.RNA.markers.down.txt",row.names=FALSE,sep="\t", quote = FALSE)


p.up=DoHeatmap(hiv, features = up.genes, group.by = 'group') + NoLegend()
p.down=DoHeatmap(hiv, features = down.genes, group.by = 'group') + NoLegend()
p.up
p.down

dev.off()

pdf("pdf/four/figure1.3.pdf", width = 5, height = 8)
DoHeatmap(hiv, features = rownames(SAHA.RNA.markers.down)[1:10], group.by = 'group') + NoLegend()
DoHeatmap(hiv, features = rownames(Pros.RNA.markers.down)[1:10], group.by = 'group') + NoLegend()
DoHeatmap(hiv, features = rownames(iBET.RNA.markers.down)[1:10], group.by = 'group') + NoLegend()

DoHeatmap(hiv, features = rownames(SAHA.RNA.markers.up)[1:10], group.by = 'group') + NoLegend()
DoHeatmap(hiv, features = rownames(Pros.RNA.markers.up)[1:10], group.by = 'group') + NoLegend()
DoHeatmap(hiv, features = rownames(iBET.RNA.markers.up)[1:10], group.by = 'group') + NoLegend()

dev.off()

intersect(a,intersect(b,c))
#"LTB"  "TC2N"


all.markers = rbind(Pros.RNA.markers,rbind(iBET.RNA.markers,SAHA.RNA.markers))
all.markers= data.frame(all.markers)
all.markers$gene = rownames(all.markers)

linked.gene <- read.csv("linkage_analysis_output_wm/linked.gene.csv")
names(linked.gene) <- c("gene", "linked.cor", "linked.pval", "linked.pval.adj")
names(all.markers) <- c("DEG.p_val", "DEG.avg_log2FC", "DEG.pct.1", "DEG.pct.2", "DEG.p_val_adj","gene")

linked.gene.vs.DEG <- merge(linked.gene, all.markers, by = "gene")
compare <- linked.gene.vs.DEG[order(linked.gene.vs.DEG[,4]),c(1,2,3,4,5,9)]
options(width = 300)
compare$log.linked.pval = sqrt(-log(compare$linked.pval.adj))
compare$log.DEG_pval = sqrt(-log(compare$DEG.p_val_adj))
compare1 <- do.call(data.frame, lapply(compare,function(x) replace(x, is.infinite(x), 30)))
compare2 = data.frame(compare1[,c(7,8)])

pdf("pdf/four/figure1.3.pdf", width = 15, height = 15)
library(ggpubr)

ggscatter(compare2, x = "log.DEG_pval", y = "log.linked.pval",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "sqrt(-log(compare$DEG.p_val_adj))", ylab = "sqrt(-log(compare$linked.pval.adj))")


dev.off()

load("motif.mat.rda")

DATF=matrix(nrow=nrow(motif.mat), ncol=7)
colnames(DATF)=c('Mean.DMSO','Mean.SAHA','Mean.Prostratin','Mean.iBET151','SAHA','Prostratin','iBET151')
DATF=as.data.frame(DATF)

for(i in 1:nrow(DATF)){
  motif.i=motif.mat[i,]
  DATF[i,1]=mean(motif.i[hiv$group=='DMSO'])
  DATF[i,2]=mean(motif.i[hiv$group=='SAHA'])
  DATF[i,3]=mean(motif.i[hiv$group=='Prostratin'])
  DATF[i,4]=mean(motif.i[hiv$group=='iBET151'])
  DATF[i,5]=wilcox.test(motif.i[hiv$group=='SAHA'],motif.i[hiv$group=='DMSO'])$p.value
  DATF[i,6]=wilcox.test(motif.i[hiv$group=='Prostratin'],motif.i[hiv$group=='DMSO'])$p.value
  DATF[i,7]=wilcox.test(motif.i[hiv$group=='iBET151'],motif.i[hiv$group=='DMSO'])$p.value
}

SAHA.adj=qvalue(DATF$SAHA)$qvalues
Prostratin.adj=qvalue(DATF$Prostratin,pi0 = 1)$qvalues
iBET151.adj=qvalue(DATF$iBET151)$qvalues


DATF=cbind(DATF, SAHA.adj, Prostratin.adj, iBET151.adj)
rownames(DATF)=rownames(motif.mat)

pdf("pdf/four/figure2.pdf", width = 10, height = 10)
x=list(SAHA.down=which((DATF$Mean.SAHA<DATF$Mean.DMSO) & (DATF$SAHA.adj<0.05)),
       Prostratin.down=which((DATF$Mean.Prostratin<DATF$Mean.DMSO) & (DATF$Prostratin.adj<0.05)),
       iBET151.down=which((DATF$Mean.iBET151<DATF$Mean.DMSO) & (DATF$iBET151.adj<0.05)))
p1=ggvenn(x,
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)
p1

x=list(SAHA.up=which((DATF$Mean.SAHA>DATF$Mean.DMSO) & (DATF$SAHA.adj<0.05)),
       Prostratin.up=which((DATF$Mean.Prostratin>DATF$Mean.DMSO) & (DATF$Prostratin.adj<0.05)),
       iBET151.up=which((DATF$Mean.iBET151>DATF$Mean.DMSO) & (DATF$iBET151.adj<0.05)))
p2=ggvenn(x,
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)
p2

p1|p2
dev.off()


DATF <-  read.csv('../HIV1/downsample/linkage_analysis_output_wm/DATF.csv')
DATF <-  read.csv('linkage_analysis_output_wm/DATF.csv')
SAHA.down=which((DATF$Mean.SAHA<DATF$Mean.DMSO) & (DATF$SAHA.pval.adj<0.05))
Prostratin.down=which((DATF$Mean.Prostratin<DATF$Mean.DMSO) & (DATF$Prostratin.pval.adj<0.05))
iBET151.down=which((DATF$Mean.iBET151<DATF$Mean.DMSO) & (DATF$iBET151.pval.adj<0.05))



intersect.down=intersect(intersect(SAHA.down, Prostratin.down), iBET151.down)


SAHA.up=which((DATF$Mean.SAHA>DATF$Mean.DMSO) & (DATF$SAHA.pval.adj<0.05))
Prostratin.up=which((DATF$Mean.Prostratin>DATF$Mean.DMSO) & (DATF$Prostratin.pval.adj<0.05))
iBET151.up=which((DATF$Mean.iBET151>DATF$Mean.DMSO) & (DATF$iBET151.pval.adj<0.05))

intersect.up=intersect(intersect(SAHA.up, Prostratin.up), iBET151.up)
colnames(DATF)=c('Mean.DMSO','Mean.SAHA',       'Mean.Prostratin','Mean.iBET151',
                 'SAHA.pval',   'Prostratin.pval',      'iBET151.pval',
                 'SAHA.pval.adj','Prostratin.pval.adj', 'iBET151.pval.adj')

cauchy.comb=function(Pvals){
  Pvals[Pvals>0.999]=0.999
  is.small<-(Pvals<1e-15)
  Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
  Pvals[is.small]<-1/Pvals[is.small]/pi
  cct.stat<-mean(Pvals)
  pval<-pcauchy(cct.stat,lower.tail = F)
  return(pval)
}
load("TF.rda")
#TF=ConvertMotifID(hiv, id=rownames(motif.mat))
DATF=cbind(TF, motif=rownames(DATF),DATF)
DATF.combined.pval=apply(DATF[,c("SAHA.pval","Prostratin.pval","iBET151.pval" )],1, cauchy.comb)
DATF.combined.pval.adj=qvalue(DATF.combined.pval,pi0 = 1)$qvalues
DATF=cbind(DATF, DATF.combined.pval, DATF.combined.pval.adj)

head(DATF)
write.csv(DATF, file='linkage_analysis_output_wm/DATF.csv', row.names = F)
DATF <- read.csv('linkage_analysis_output_wm/DATF.csv',header = T)


new_assay <- CreateAssayObject(counts = motif.mat)
hiv[["chromvar"]] <- new_assay
DefaultAssay(hiv) <- "chromvar"

pdf("pdf/four/figure3.pdf", width = 10, height = 10)

par(mfrow=c(2,2))
hist(DATF$SAHA.pval, xlab='p-value', main='DATF SAHA')
hist(DATF$Prostratin.pval, xlab='p-value', main='DATF Prostratin')
hist(DATF$iBET151.pval, xlab='p-value', main='DATF iBET151')
hist(DATF$DATF.combined.pval, xlab='p-value', main='DATF Cauchy Combination')

dev.off()


pdf("pdf/four/figure4.1.pdf", width = 15, height = 13)
data=data.frame(motif_deviation=hiv@assays$chromvar@data[intersect(intersect.down2,intersect.down1)[10],],group=hiv$group)
p1 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[intersect(intersect.down2,intersect.down1)[10]],': combined adj. p-val =', signif(DATF[intersect(intersect.down2,intersect.down1)[10], 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")
data=data.frame(motif_deviation=hiv@assays$chromvar@data[intersect(intersect.down2,intersect.down1)[11],],group=hiv$group)
p2 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[intersect(intersect.down2,intersect.down1)[11]],': combined adj. p-val =', signif(DATF[intersect(intersect.down2,intersect.down1)[11], 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")

data=data.frame(motif_deviation=hiv@assays$chromvar@data[intersect(intersect.down2,intersect.down1)[7],],group=hiv$group)
p3 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[intersect(intersect.down2,intersect.down1)[7]],': combined adj. p-val =', signif(DATF[intersect(intersect.down2,intersect.down1)[7], 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")

data=data.frame(motif_deviation=hiv@assays$chromvar@data[intersect(intersect.down2,intersect.down1)[4],],group=hiv$group)
p4 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[intersect(intersect.down2,intersect.down1)[4]],': combined adj. p-val =', signif(DATF[intersect(intersect.down2,intersect.down1)[4], 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")

data=data.frame(motif_deviation=hiv@assays$chromvar@data[intersect(intersect.down2,intersect.down1)[5],],group=hiv$group)
p5 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[intersect(intersect.down2,intersect.down1)[5]],': combined adj. p-val =', signif(DATF[intersect(intersect.down2,intersect.down1)[5], 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")

data=data.frame(motif_deviation=hiv@assays$chromvar@data[intersect(intersect.down2,intersect.down1)[6],],group=hiv$group)
p6 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[intersect(intersect.down2,intersect.down1)[6]],': combined adj. p-val =', signif(DATF[intersect(intersect.down2,intersect.down1)[6], 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")

p1 + p2 + p3 +p4 +p5+p6 + plot_layout(ncol=3)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("pdf/four/figure5.1.pdf", width = 15, height = 13)
data=data.frame(motif_deviation=hiv@assays$chromvar@data[intersect(intersect.up2,intersect.up1)[1],],group=hiv$group)
p1 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[intersect(intersect.up2,intersect.up1)[1]],': combined adj. p-val =', signif(DATF[intersect(intersect.up2,intersect.up1)[1], 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")
data=data.frame(motif_deviation=hiv@assays$chromvar@data[intersect(intersect.up2,intersect.up1)[2],],group=hiv$group)
p2 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[intersect(intersect.up2,intersect.up1)[2]],': combined adj. p-val =', signif(DATF[intersect(intersect.up2,intersect.up1)[2], 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")

data=data.frame(motif_deviation=hiv@assays$chromvar@data[intersect(intersect.up2,intersect.up1)[3],],group=hiv$group)
p3 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[intersect(intersect.up2,intersect.up1)[3]],': combined adj. p-val =', signif(DATF[intersect(intersect.up2,intersect.up1)[3], 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")

data=data.frame(motif_deviation=hiv@assays$chromvar@data[intersect(intersect.up2,intersect.up1)[4],],group=hiv$group)
p4 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[intersect(intersect.up2,intersect.up1)[4]],': combined adj. p-val =', signif(DATF[intersect(intersect.up2,intersect.up1)[4], 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")

data=data.frame(motif_deviation=hiv@assays$chromvar@data[intersect(intersect.up2,intersect.up1)[5],],group=hiv$group)
p5 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[intersect(intersect.up2,intersect.up1)[5]],': combined adj. p-val =', signif(DATF[intersect(intersect.up2,intersect.up1)[5], 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")

data=data.frame(motif_deviation=hiv@assays$chromvar@data[intersect(intersect.up2,intersect.up1)[6],],group=hiv$group)
p6 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[intersect(intersect.up2,intersect.up1)[6]],': combined adj. p-val =', signif(DATF[intersect(intersect.up2,intersect.up1)[6], 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")

p1 + p2 + p3 +p4 +p5+p6 + plot_layout(ncol=3)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("pdf/four/figure4.pdf", width = 15, height = 13)
#DONOR 2, 24, 27, 5 
p1=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.down2,intersect.down1)[10]], pt.size = 0)+
  ggtitle(paste(TF[intersect(intersect.down,intersect.down1)[10]],': combined adj. p-val =', signif(DATF[intersect(intersect.down,intersect.down1)[10], 'DATF.combined.pval.adj'],3)))
p2=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.down2,intersect.down1)[11]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.down,intersect.down1)[11]],': combined adj. p-val =', signif(DATF[intersect(intersect.down,intersect.down1)[11], 'DATF.combined.pval.adj'],3)))
p3=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.down2,intersect.down1)[7]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.down,intersect.down1)[7]],': combined adj. p-val =', signif(DATF[intersect(intersect.down,intersect.down1)[7], 'DATF.combined.pval.adj'],3)))
p4=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.down2,intersect.down1)[4]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.down,intersect.down1)[4]],': combined adj. p-val =', signif(DATF[intersect(intersect.down,intersect.down1)[4], 'DATF.combined.pval.adj'],3)))
p5=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.down2,intersect.down1)[5]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.down,intersect.down1)[5]],': combined adj. p-val =', signif(DATF[intersect(intersect.down,intersect.down1)[4], 'DATF.combined.pval.adj'],3)))
p6=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.down2,intersect.down1)[6]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.down,intersect.down1)[6]],': combined adj. p-val =', signif(DATF[intersect(intersect.down,intersect.down1)[4], 'DATF.combined.pval.adj'],3)))
#p7=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[7]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.down,intersect.down1)[7]],': combined adj. p-val =', signif(DATF[intersect.down[4], 'DATF.combined.pval.adj'],3)))
#p8=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[8]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.down,intersect.down1)[8]],': combined adj. p-val =', signif(DATF[intersect.down[4], 'DATF.combined.pval.adj'],3)))
#p9=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[9]], pt.size = 0)+ggtitle(paste(TF[intersect.down[9]],': combined adj. p-val =', signif(DATF[intersect.down[4], 'DATF.combined.pval.adj'],3)))
#p10=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[10]], pt.size = 0)+ggtitle(paste(TF[intersect.down[10]],': combined adj. p-val =', signif(DATF[intersect.down[4], 'DATF.combined.pval.adj'],3)))
#p11=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[11]], pt.size = 0)+ggtitle(paste(TF[intersect.down[11]],': combined adj. p-val =', signif(DATF[intersect.down[4], 'DATF.combined.pval.adj'],3)))
#p12=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[12]], pt.size = 0)+ggtitle(paste(TF[intersect.down[12]],': combined adj. p-val =', signif(DATF[intersect.down[4], 'DATF.combined.pval.adj'],3)))
#p13=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[13]])+ggtitle(paste(TF[intersect.down[13]],': combined adj. p-val =', signif(DATF[intersect.down[4], 'DATF.combined.pval.adj'],3)))
#p14=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[14]])+ggtitle(paste(TF[intersect.down[14]],': combined adj. p-val =', signif(DATF[intersect.down[4], 'DATF.combined.pval.adj'],3)))
#p15=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[15]])+ggtitle(paste(TF[intersect.down[15]],': combined adj. p-val =', signif(DATF[intersect.down[4], 'DATF.combined.pval.adj'],3)))
#p16=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[16]])+ggtitle(paste(TF[intersect.down[16]],': combined adj. p-val =', signif(DATF[intersect.down[4], 'DATF.combined.pval.adj'],3)))

p1 + p2 + p3 +p4 +p5+p6 + plot_layout(ncol=3)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))

dev.off()

pdf("pdf/four/figure5.pdf", width = 15, height = 13)
p1=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.up,intersect.up1)[1]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.up,intersect.up1)[1]],': combined adj. p-val =', signif(DATF[intersect(intersect.up,intersect.up1)[1], 'DATF.combined.pval.adj'],3)))
p2=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.up,intersect.up1)[2]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.up,intersect.up1)[2]],': combined adj. p-val =', signif(DATF[intersect(intersect.up,intersect.up1)[2], 'DATF.combined.pval.adj'],3)))
p3=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.up,intersect.up1)[3]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.up,intersect.up1)[3]],': combined adj. p-val =', signif(DATF[intersect(intersect.up,intersect.up1)[3], 'DATF.combined.pval.adj'],3)))
p4=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.up,intersect.up1)[4]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.up,intersect.up1)[4]],': combined adj. p-val =', signif(DATF[intersect(intersect.up,intersect.up1)[4], 'DATF.combined.pval.adj'],3)))
p5=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.up,intersect.up1)[5]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.up,intersect.up1)[5]],': combined adj. p-val =', signif(DATF[intersect(intersect.up,intersect.up1)[5], 'DATF.combined.pval.adj'],3)))
p6=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.up,intersect.up1)[6]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.up,intersect.up1)[6]],': combined adj. p-val =', signif(DATF[intersect(intersect.up,intersect.up1)[6], 'DATF.combined.pval.adj'],3)))
#p7=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.up,intersect.up1)[7]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.up,intersect.up1)[7]],': combined adj. p-val =', signif(DATF[intersect(intersect.up,intersect.up1)[7], 'DATF.combined.pval.adj'],3)))
#p8=VlnPlot(hiv, features = rownames(motif.mat)[intersect(intersect.up,intersect.up1)[8]], pt.size = 0)+ggtitle(paste(TF[intersect(intersect.up,intersect.up1)[8]],': combined adj. p-val =', signif(DATF[intersect(intersect.up,intersect.up1)[8], 'DATF.combined.pval.adj'],3)))
#p9=VlnPlot(hiv, features = rownames(motif.mat)[intersect.up[9]], pt.size = 0)+ggtitle(paste(TF[intersect.up[9]],': combined adj. p-val =', signif(DATF[intersect.up[4], 'DATF.combined.pval.adj'],3)))
#p10=VlnPlot(hiv, features = rownames(motif.mat)[intersect.up[10]], pt.size = 0)+ggtitle(paste(TF[intersect.up[10]],': combined adj. p-val =', signif(DATF[intersect.up[4], 'DATF.combined.pval.adj'],3)))
#p11=VlnPlot(hiv, features = rownames(motif.mat)[intersect.up[11]], pt.size = 0)+ggtitle(paste(TF[intersect.up[11]],': combined adj. p-val =', signif(DATF[intersect.up[4], 'DATF.combined.pval.adj'],3)))
#p12=VlnPlot(hiv, features = rownames(motif.mat)[intersect.up[12]], pt.size = 0)+ggtitle(paste(TF[intersect.up[12]],': combined adj. p-val =', signif(DATF[intersect.up[4], 'DATF.combined.pval.adj'],3)))
#p13=VlnPlot(hiv, features = rownames(motif.mat)[intersect.up[13]])+ggtitle(paste(TF[intersect.up[13]],': combined adj. p-val =', signif(DATF[intersect.up[4], 'DATF.combined.pval.adj'],3)))
#p14=VlnPlot(hiv, features = rownames(motif.mat)[intersect.up[14]])+ggtitle(paste(TF[intersect.up[14]],': combined adj. p-val =', signif(DATF[intersect.up[4], 'DATF.combined.pval.adj'],3)))
#p15=VlnPlot(hiv, features = rownames(motif.mat)[intersect.up[15]])+ggtitle(paste(TF[intersect.up[15]],': combined adj. p-val =', signif(DATF[intersect.up[4], 'DATF.combined.pval.adj'],3)))
#p16=VlnPlot(hiv, features = rownames(motif.mat)[intersect.up[16]])+ggtitle(paste(TF[intersect.up[16]],': combined adj. p-val =', signif(DATF[intersect.up[4], 'DATF.combined.pval.adj'],3)))

#p1 + p2 + plot_layout(ncol=2)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p1 + p2 + p3 +p4 +p5+p6 + plot_layout(ncol=3)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))


dev.off()

pdf("pdf/four/figure6.pdf", width = 15, height = 13)
data=data.frame(motif_deviation=hiv@assays$chromvar@data[which(TF == "JUN(var.2)"),],group=hiv$group)
p1 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[which(TF == "JUN")],': combined adj. p-val =', signif(DATF2[which(TF == "JUN"), 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")

data=data.frame(motif_deviation=hiv@assays$chromvar@data[which(TF == "FOS"),],group=hiv$group)
p2 <- ggplot(data, aes(x=group, y=motif_deviation, fill=group)) + geom_boxplot()  + ggtitle(paste(TF[which(TF == "FOS")],': combined adj. p-val =', signif(DATF2[which(TF == "FOS"), 'DATF.combined.pval.adj'],3))) + theme(legend.position="none")

p1|p2
dev.off()



# DATF
head(DATF)
colnames(DATF)[3:(ncol(DATF)-2)]=paste0('DATF.', colnames(DATF)[3:(ncol(DATF)-2)])

# Linked TF motif accessibility
linkedTF=read.csv('linkage_analysis_output_wm/linked.TF.csv')
#dim = 633
linkedTF=linkedTF[match(DATF$TF,linkedTF$TF),]
colnames(linkedTF)[3:ncol(linkedTF)]=paste0('linkedTF.',colnames(linkedTF)[3:ncol(linkedTF)])
linked.TF <- linkedTF
#dim = 633


# Linked TF expression
linkedTF.exp.temp=as.data.frame(read.csv('linkage_analysis_output_wm/linked.TF.gene.csv'))
#dim = 176 4
temp <- read.delim("linkage_analysis_output_wm/linked.TF.gene.output.txt",header = TRUE, sep = "\t")


linkedTF.exp=matrix(nrow=nrow(DATF), ncol=ncol(linkedTF.exp.temp))
colnames(linkedTF.exp)=colnames(linkedTF.exp.temp)
rownames(linkedTF.exp)=rownames(DATF)
linkedTF.exp=as.data.frame(linkedTF.exp)

linkedTF.exp[match(linkedTF.exp.temp$gene,DATF$TF),]=linkedTF.exp.temp
rm(linkedTF.exp.temp)
colnames(linkedTF.exp)[2:ncol(linkedTF.exp)]=paste0('linkedTF.exp.',colnames(linkedTF.exp)[2:ncol(linkedTF.exp)])


#TF with enriched motif in linked peaks

enrichedTF=read.csv('linkage_analysis_output_wm/enriched.motif.csv')
enrichedTF=enrichedTF[match(DATF$TF,enrichedTF$motif.name),]
enrichedTF = enrichedTF[,c(1:7,9,8)]
colnames(enrichedTF)[2:(ncol(enrichedTF)-1)]=paste0('enrichedTF.',colnames(enrichedTF)[2:(ncol(enrichedTF)-1)])

TF.output=cbind(DATF, linkedTF, enrichedTF, linkedTF.exp)
TF.output=TF.output[,-c(15,16,20,27,29)]
write.csv(TF.output, file='linkage_analysis_output_wm/TF.output.csv', row.names = F)

linkedTF.exp = cbind(linkedTF.exp,DATF[,2])

pdf("pdf/five/figure1.pdf", width = 10, height = 10)
library(ggvenn)
x=list(linkedTF=linked.TF[linked.TF$linkedTF.pval.adj<=0.1,2],
       enrichedTF=enrichedTF[enrichedTF$enrichedTF.pvalue<=0.1,1],
       DATF=DATF[DATF$DATF.combined.pval.adj<0.1,2])
p1=ggvenn(x,
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)
p1




x=list(linkedTF=linked.TF[linked.TF$linkedTF.pval.adj<=0.1,2],
       enrichedTF=enrichedTF[enrichedTF$enrichedTF.pvalue<=0.1,1],
       DATF=DATF[DATF$DATF.combined.pval.adj<0.1,2],
       linked.TF.exp=linkedTF.exp[linkedTF.exp$linkedTF.exp.pval.adj<=0.1,5])
p2=ggvenn(x,
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)
p2
dev.off()


overlapped.motifs=intersect(intersect(linked.TF[linked.TF$linkedTF.pval.adj<=0.1,2], enrichedTF[enrichedTF$enrichedTF.pvalue<=0.1,1]), DATF[DATF$DATF.combined.pval.adj<0.1,2])
overlapped.motifs
linked.TF[match(overlapped.motifs, linked.TF[,2]),2]

head(linked.TF[match(overlapped.motifs, linked.TF[,2]),],20)
enrichedTF[match(head(linked.TF[match(overlapped.motifs, linked.TF[,2]),],20)$TF, enrichedTF$motif.name),-1]
DATF[match(head(linked.TF[match(overlapped.motifs, linked.TF[,2]),],20)$TF, DATF$TF),]

# Below is to remove FOS and JUN
enrichedTF.filter=enrichedTF[!grepl('JUN', enrichedTF$motif.name) & !grepl('FOS', enrichedTF$motif.name) &
                               !grepl('BATF', enrichedTF$motif.name) & !grepl('BACH', enrichedTF$motif.name) &
                               !grepl('JDP2', enrichedTF$motif.name) & !grepl('NFE', enrichedTF$motif.name),]

overlapped.motifs=intersect(intersect(linked.TF[linked.TF$linkedTF.pval.adj<=0.1,2], enrichedTF.filter[enrichedTF.filter$enrichedTF.pvalue<=0.1,1]), DATF[DATF$DATF.combined.pval.adj<0.1,2])
overlapped.motifs
head(linked.TF[match(overlapped.motifs, linked.TF[,2]),],20)
enrichedTF[match(head(linked.TF[match(overlapped.motifs, linked.TF[,2]),],20)$TF, enrichedTF$motif.name),-1]
DATF[match(head(linked.TF[match(overlapped.motifs, linked.TF[,2]),],20)$TF, DATF$TF),]
#enrichedTF.filter[match(overlapped.motifs, enrichedTF[,1]),]
#DATF[match(overlapped.motifs, DATF[,2]),]





AP1.genes <- c("BATF3", "BATF2", "BATF", "JDP2","ATF7", "ATF6", "ATF5", "ATF4", "ATF3", "ATF2", "ATF1", "JUND",
               "JUNB","cJUN", "FOSL2", "FOSL1","FOSB", "cFOS")
result <- filter(TF.output, grepl(paste(AP1.genes, collapse="|"), TF))
result[,c(1,2,15:17)]
result[,c(1:14)]
result[,c(1,2,18:23)]
na.omit(result[,c(1,2,24:26)])


load("motif.mat.rda")

load("pRNA.rda")

pdf("pdf/five/figure1.5.pdf", width = 15, height = 15)
# Plot for top linked TFs
par(mfrow=c(3,4))
for(TF.name in AP1.genes[AP1.genes %in% TF]){
  temp=summary(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]))
  smoothScatter(motif.mat[which(TF==TF.name),],sqrt(pRNA), xlab='TF motif score.', ylab='sqrt(percentage of HIV RNA reads)',
                main=paste(TF.name, '\np-val =',signif(temp$coefficients[2,4],3), 'r =', signif(cor(sqrt(pRNA),motif.mat[which(TF==TF.name),]),3)))
  abline(lm(sqrt(pRNA)~motif.mat[which(TF==TF.name),]), lty=2, col='black')
}
par(mfrow=c(1,1))
dev.off()





####################################
#### TF network analysis
####################################

# Release the treshold min_cells, argument passed down to vst in sctransform
hiv <- SCTransform(hiv, verbose = FALSE, return.only.var.genes = FALSE, min_cells=1) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

TF.output <- read.csv("linkage_analysis_output_wm/TF.output.csv",header = TRUE, sep = ",")

# Focus on TFs called by at least two testing
temp=TF.output$linkedTF.exp.pval.adj
temp[is.na(temp)]=1
alpha=0.01
TF.sel=TF.output$TF[(as.numeric(TF.output$DATF.combined.pval.adj<=alpha)+
                       as.numeric(TF.output$linkedTF.pval.adj<=alpha)+
                       as.numeric(TF.output$enrichedTF.pvalue<=alpha)+
                       as.numeric(temp<=alpha))>=2]
length(TF.sel)
TF.sel[grep('SNA', TF.sel)]

TF.sel=intersect(TF.sel, rownames(hiv@assays$SCT))

# TF exp
TF.exp=hiv@assays$SCT@data[match(TF.sel, rownames(hiv@assays$SCT)),]
# TF motif score
TF.access=hiv@assays$chromvar@data[match(TF.sel, TF),]
length(TF.sel)
dim(TF.exp)
dim(TF.access)

TF.exp.pc=TF.exp%*%svd(TF.exp)$v[,1:10]
TF.access.pc=TF.access%*%svd(TF.access)$v[,1:10]

TF.exp.access.pc=cbind(TF.exp.pc, TF.access.pc)
TF.consensus.pc=as.matrix(TF.exp.access.pc%*%svd(TF.exp.access.pc)$v)
TF.cor=cor(t(TF.consensus.pc))


library(pheatmap)
pdf("pdf/five/figure2.pdf", width = 20, height = 20)
pheatmap(TF.cor)
dev.off()



library(caret)
library(Rtsne)
set.seed(1)

tsne_model_1 = Rtsne(as.matrix(TF.cor), check_duplicates=FALSE, pca=F, perplexity=30, theta=0.5, dims=2)
d_tsne_1 = as.data.frame(tsne_model_1$Y)

tsne_model_2 = Rtsne(as.matrix(TF.cor), check_duplicates=FALSE, pca=T, perplexity=30, theta=0.5, dims=2)
d_tsne_2 = as.data.frame(tsne_model_2$Y)

pdf("hiv/pdf/seven/figure2.pdf", width = 10, height = 10)
ggplot(d_tsne_2, aes(x=V1, y=V2)) +
  geom_point(size=1) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("X1") + ylab("X2") +
  ggtitle("t-SNE & PCA") +
  theme_light(base_size=20)
dev.off()

pdf("hiv/pdf/seven/figure2.1.pdf", width = 10, height = 10)
ggplot(d_tsne_1, aes(x=V1, y=V2)) +
  geom_point(size=1) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("X1") + ylab("X2") +
  ggtitle("t-SNE") +
  theme_light(base_size=20)
dev.off()



RunPCA(d_tsne_1)

library(PCAtools)



library(psych)

pc <- prcomp(d_tsne_1,
             center = TRUE,
             scale. = TRUE)
attributes(pc)








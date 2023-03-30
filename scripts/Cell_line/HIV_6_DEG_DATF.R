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

load('processed_data/hiv_predict.rda')

# Differentially expressed RNA between conditions

DefaultAssay(hiv)='SCT'
table(hiv$group)

SAHA.RNA.markers <- FindMarkers(hiv, ident.1 = 'SAHA', ident.2='DMSO',min.pct = 0.25)
SAHA.RNA.markers=SAHA.RNA.markers[SAHA.RNA.markers$p_val_adj<=0.05,]
SAHA.RNA.markers.down = SAHA.RNA.markers[SAHA.RNA.markers$avg_log2FC<0,]
SAHA.RNA.markers.up = SAHA.RNA.markers[SAHA.RNA.markers$avg_log2FC>0,]

write.csv(SAHA.RNA.markers.down, file='output/supp_table1_SAHA.RNA.markers.down.csv', row.names = T, quote = F)
write.csv(SAHA.RNA.markers.up, file='output/supp_table1_SAHA.RNA.markers.up.csv', row.names = T, quote = F)

Pros.RNA.markers <- FindMarkers(hiv, ident.1 = 'Prostratin', ident.2='DMSO',min.pct = 0.25)
Pros.RNA.markers=Pros.RNA.markers[Pros.RNA.markers$p_val_adj<=0.05,]
Pros.RNA.markers.down = Pros.RNA.markers[Pros.RNA.markers$avg_log2FC<0,]
Pros.RNA.markers.up = Pros.RNA.markers[Pros.RNA.markers$avg_log2FC>0,]

write.csv(Pros.RNA.markers.down, file='output/supp_table1_Pros.RNA.markers.down.csv', row.names = T, quote = F)
write.csv(Pros.RNA.markers.up, file='output/supp_table1_Pros.RNA.markers.up.csv', row.names = T, quote = F)

iBET.RNA.markers <- FindMarkers(hiv, ident.1 = 'iBET151', ident.2='DMSO',min.pct = 0.25)
iBET.RNA.markers=iBET.RNA.markers[iBET.RNA.markers$p_val_adj<=0.05,]
iBET.RNA.markers.down = iBET.RNA.markers[iBET.RNA.markers$avg_log2FC<0,]
iBET.RNA.markers.up = iBET.RNA.markers[iBET.RNA.markers$avg_log2FC>0,]

write.csv(iBET.RNA.markers.down, file='output/supp_table1_iBET.RNA.markers.down.csv', row.names = T, quote = F)
write.csv(iBET.RNA.markers.up, file='output/supp_table1_iBET.RNA.markers.up.csv', row.names = T, quote = F)

library(ggvenn)
x=list(SAHA.down=rownames(SAHA.RNA.markers.down),
       Prostratin.down=rownames(Pros.RNA.markers.down),
       iBET151.down=rownames(iBET.RNA.markers.down))
p1=ggvenn(x, 
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)
p1

x=list(SAHA.up=rownames(SAHA.RNA.markers.up),
       Prostratin.up=rownames(Pros.RNA.markers.up),
       iBET151.up=rownames(iBET.RNA.markers.up))
p2=ggvenn(x, 
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)
p2
intersect(intersect(rownames(SAHA.RNA.markers.up), rownames(Pros.RNA.markers.up)), 
          rownames(iBET.RNA.markers.up))

p=p1+p2+plot_layout(ncol=1)
p
ggsave(p, file='output/figure3_venn.pdf', width=3, height=6)

up.genes=unique(c(rownames(SAHA.RNA.markers.up)[1:20],
                  rownames(Pros.RNA.markers.up)[1:20],
                  rownames(iBET.RNA.markers.up)[1:20]))
up.genes=up.genes[1:50]
p.up=DoHeatmap(hiv, features = up.genes, group.by = 'group')
ggsave('output/figure3_heatmap.up.legend.pdf', plot=p.up, width = 5, height=8)

down.genes=unique(c(rownames(SAHA.RNA.markers.down)[1:50],
                  rownames(Pros.RNA.markers.down)[1:50],
                  rownames(iBET.RNA.markers.down)[1:50]))
down.genes=down.genes[1:50]
p.down=DoHeatmap(hiv, features = down.genes, group.by = 'group')
ggsave('output/figure3_heatmap.down.legend.pdf', plot=p.down, width = 5, height=8)



# Differentially accessible TFs
load('processed_data/hiv_chromvar.rda')

# # Below is testing differentially accessible peaks
# # This does not make too much sense because the data is too sparse
# DefaultAssay(hiv)='ATAC'
# SAHA.da_peaks <- FindMarkers(
#   object = hiv,
#   ident.1 = 'SAHA', ident.2='DMSO',
#   test.use = 'LR',
#   min.pct = 0.005,
#   logfc.threshold = 0.1
# )
# 
# Pros.da_peaks <- FindMarkers(
#   object = hiv,
#   ident.1 = 'Prostratin', ident.2='DMSO',
#   test.use = 'LR',
#   min.pct = 0.005,
#   logfc.threshold = 0.1
# )
# 
# iBET.da_peaks <- FindMarkers(
#   object = hiv,
#   ident.1 = 'iBET151', ident.2='DMSO',
#   test.use = 'LR',
#   min.pct = 0.005,
#   logfc.threshold = 0.1
# )

# Differentially accessible TFs
dim(hiv@assays$chromvar@data)
dim(motif.mat)
length(hiv$group)

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
Prostratin.adj=qvalue(DATF$Prostratin)$qvalues
iBET151.adj=qvalue(DATF$iBET151)$qvalues


DATF=cbind(DATF, SAHA.adj, Prostratin.adj, iBET151.adj)
rownames(DATF)=rownames(motif.mat)


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

p=p1|p2
ggsave(p, file='output/figure4_DATF_venn.pdf', width=6, height=3)


SAHA.down=which((DATF$Mean.SAHA<DATF$Mean.DMSO) & (DATF$SAHA.adj<0.05))
Prostratin.down=which((DATF$Mean.Prostratin<DATF$Mean.DMSO) & (DATF$Prostratin.adj<0.05))
iBET151.down=which((DATF$Mean.iBET151<DATF$Mean.DMSO) & (DATF$iBET151.adj<0.05))

intersect.down=intersect(intersect(SAHA.down, Prostratin.down), iBET151.down)

# library(pheatmap)
# pheatmap(pmin(2,pmax(-2,motif.mat[down.TF,order(hiv$group)])), show_colnames = F, cluster_rows = F, cluster_cols = F)

colnames(DATF)=c('Mean.DMSO','Mean.SAHA',	'Mean.Prostratin','Mean.iBET151',
                 'SAHA.pval',	'Prostratin.pval',	'iBET151.pval',
                 'SAHA.pval.adj','Prostratin.pval.adj',	'iBET151.pval.adj')

cauchy.comb=function(Pvals){
  Pvals[Pvals>0.999]=0.999
  is.small<-(Pvals<1e-15)
  Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
  Pvals[is.small]<-1/Pvals[is.small]/pi
  cct.stat<-mean(Pvals)
  pval<-pcauchy(cct.stat,lower.tail = F)
  return(pval)
}

DATF=cbind(TF, motif=rownames(DATF),DATF)
DATF.combined.pval=apply(DATF[,c("SAHA.pval","Prostratin.pval","iBET151.pval" )],1, cauchy.comb)
DATF.combined.pval.adj=qvalue(DATF.combined.pval)$qvalues
DATF=cbind(DATF, DATF.combined.pval, DATF.combined.pval.adj)

# Look at distributions of nominal p-values
par(mfrow=c(2,2))
hist(DATF$SAHA.pval, xlab='p-value', main='DATF SAHA')
hist(DATF$Prostratin.pval, xlab='p-value', main='DATF Prostratin')
hist(DATF$iBET151.pval, xlab='p-value', main='DATF iBET151')
hist(DATF$DATF.combined.pval, xlab='p-value', main='DATF Cauchy Combination')


p1=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[1]])+ggtitle(paste(TF[intersect.down[1]],': p.cauchy.adj =', signif(DATF[intersect.down[1], 'DATF.combined.pval.adj'],3)))
p2=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[2]])+ggtitle(paste(TF[intersect.down[2]],': p.cauchy.adj =', signif(DATF[intersect.down[2], 'DATF.combined.pval.adj'],3)))
p3=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[3]])+ggtitle(paste(TF[intersect.down[3]],': p.cauchy.adj =', signif(DATF[intersect.down[3], 'DATF.combined.pval.adj'],3)))
p4=VlnPlot(hiv, features = rownames(motif.mat)[intersect.down[4]])+ggtitle(paste(TF[intersect.down[4]],': p.cauchy.adj =', signif(DATF[intersect.down[4], 'DATF.combined.pval.adj'],3)))

p = p1 + p2 + p3 +p4 +plot_layout(ncol=2)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p
ggsave(p, file='output/figure4_four_intersected_TFs.pdf', width=8, height=8)

head(DATF)
write.csv(DATF, file='linkage_analysis_output/DATF.csv', row.names = F)

# Let's look at AP-1
# All TFs that contain either "JUN" or "FOS"
TF[unique(c(grep('JUN', TF), grep('FOS', TF)))]


TF.name='FOS'
p1=VlnPlot(hiv, features = rownames(motif.mat)[which(TF==TF.name)])+ggtitle(paste(TF.name,': p.cauchy.adj =', signif(DATF[which(TF==TF.name), 'DATF.combined.pval.adj'],3)))
TF.name='FOSL1'
p2=VlnPlot(hiv, features = rownames(motif.mat)[which(TF==TF.name)])+ggtitle(paste(TF.name,': p.cauchy.adj =', signif(DATF[which(TF==TF.name), 'DATF.combined.pval.adj'],3)))
TF.name='FOSL2'
p3=VlnPlot(hiv, features = rownames(motif.mat)[which(TF==TF.name)])+ggtitle(paste(TF.name,': p.cauchy.adj =', signif(DATF[which(TF==TF.name), 'DATF.combined.pval.adj'],3)))
TF.name='JUN'
p4=VlnPlot(hiv, features = rownames(motif.mat)[which(TF==TF.name)])+ggtitle(paste(TF.name,': p.cauchy.adj =', signif(DATF[which(TF==TF.name), 'DATF.combined.pval.adj'],3)))
TF.name='JUN(var.2)'
p5=VlnPlot(hiv, features = rownames(motif.mat)[which(TF==TF.name)])+ggtitle(paste(TF.name,': p.cauchy.adj =', signif(DATF[which(TF==TF.name), 'DATF.combined.pval.adj'],3)))
TF.name='JUNB'
p6=VlnPlot(hiv, features = rownames(motif.mat)[which(TF==TF.name)])+ggtitle(paste(TF.name,': p.cauchy.adj =', signif(DATF[which(TF==TF.name), 'DATF.combined.pval.adj'],3)))
TF.name='FOS::JUN'
p7=VlnPlot(hiv, features = rownames(motif.mat)[which(TF==TF.name)])+ggtitle(paste(TF.name,': p.cauchy.adj =', signif(DATF[which(TF==TF.name), 'DATF.combined.pval.adj'],3)))
TF.name='FOS::JUN(var.2)'
p8=VlnPlot(hiv, features = rownames(motif.mat)[which(TF==TF.name)])+ggtitle(paste(TF.name,': p.cauchy.adj =', signif(DATF[which(TF==TF.name), 'DATF.combined.pval.adj'],3)))
TF.name="FOS::JUNB"
p9=VlnPlot(hiv, features = rownames(motif.mat)[which(TF==TF.name)])+ggtitle(paste(TF.name,': p.cauchy.adj =', signif(DATF[which(TF==TF.name), 'DATF.combined.pval.adj'],3)))

p1 + p2 + p3 +p4 +p5+
  p6+p7+p8+p9+plot_layout(ncol=3)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))

p=p1+p5+plot_layout(ncol=2)& NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p
ggsave(p, file='output/figure4_jun_fos.pdf', width=8, height=4)

save.image(file='processed_data/hiv_predict_DE.rda')

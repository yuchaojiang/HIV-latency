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

DATF <- read.csv("linkage_analysis_output_wm/DATF.csv")


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
temp <- read.csv("linkage_analysis_output_wm/linked.TF.gene.csv",header = TRUE)


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
colnames(enrichedTF)[2:(ncol(enrichedTF)-1)]=paste0('enrichedTF.',colnames(enrichedTF)[2:(ncol(enrichedTF)-1)])

TF.output=cbind(DATF, linkedTF, enrichedTF, linkedTF.exp)
TF.output=TF.output[,-c(15,16,20,27,28)]
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


#Compare DATF and linkedTF 

com.TF <- TF.output[,c(1,2,14,17)]
com.TF$linkedTF.transformed = sqrt(-log(com.TF$linkedTF.pval.adj))
com.TF$DATF.transformed = sqrt(-log(com.TF$DATF.combined.pval.adj))

com.TF1 <- do.call(data.frame, lapply(com.TF,function(x) replace(x, is.infinite(x), 30)))


pdf("full/pdf/five/figure1.1.pdf", width = 5, height = 5)
library(ggpubr)

ggscatter(com.TF1, x = "linkedTF.transformed", y = "DATF.transformed",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "sqrt(-log(Linked.TF.pval.adj))", ylab = "sqrt(-log(DATF.combined.pval.adj))", ylim = c(0, 30), xlim = c(0, 30))


dev.off()



AP1.genes <- c("BATF3", "BATF2", "BATF", "JDP2","ATF7", "ATF6", "ATF5", "ATF4", "ATF3", "ATF2", "ATF1", "JUND",
               "JUNB","cJUN", "FOSL2", "FOSL1","FOSB", "cFOS")
result <- filter(TF.output, grepl(paste(AP1.genes, collapse="|"), TF))
result[,c(1,2,15:17)]
result[,c(1:14)]
result[,c(1,2,18:23)]
na.omit(result[,c(1,2,24:26)])


load("full/motif.mat.rda")

load("full/pRNA.rda")

pdf("full/pdf/five/figure1.5.pdf", width = 15, height = 15)
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
  #DONOR1 alpha=0.0000000000000000000000000000001
  alpha=0.000001
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
  TF.access=motif.mat[match(TF.sel, TF),]
  length(TF.sel)
  dim(TF.exp)
  dim(TF.access)
  
  TF.exp.pc=TF.exp%*%svd(TF.exp)$v[,1:10]
  TF.access.pc=TF.access%*%svd(TF.access)$v[,1:10]
  
  TF.exp.access.pc=cbind(TF.exp.pc, TF.access.pc)
  TF.consensus.pc=as.matrix(TF.exp.access.pc%*%svd(TF.exp.access.pc)$v)
  TF.cor=cor(t(TF.consensus.pc))
  

  
  pdf("pdf/five/figure2.new.pdf", width = 20, height = 20)
  library(pheatmap)
  pheatmap(TF.cor)
  dev.off()
  
  
  TF.cor[upper.tri(TF.cor)] <- 42
  
  library(reshape2)
  TF.cor_df <- melt(TF.cor)
  library(dplyr)
  TF.cor_df <- filter(TF.cor_df, value != 42) %>% filter(Var1 != Var2)
  dim(TF.cor_df)
  
  summary(TF.cor_df$value)
  
  # create adjacency list with correlations > 0.5
  TF.adj_list <- TF.cor_df %>% filter(value > 0.5)
  names(TF.adj_list) <- c('from', 'to', 'weight')
  dim(TF.adj_list)
  
  library(igraph)
  # create igraph S3 object
  net <- graph.data.frame(TF.adj_list, directed = FALSE)
  
  # store original margins
  orig_mar <- par()$mar
  
  # set new margins to limit whitespace in plot
  par(mar=rep(.1, 4))
  
  # not much difference in the edge width given the values
  # but I included it for reference's sake
  set.seed(111)
  plot(net, layout = layout_components(net), edge.width = E(net)$weight)
  pdf("pdf/five/figure3.pdf", width = 5, height = 5)
  plot(net, layout = layout_components(net), edge.width = E(net)$weight, vertex.shape="none")
  dev.off()
  
  DATF=DATF[match(linked.TF$TF, DATF$TF),]
  enriched=enriched[match(linked.TF$TF, enriched$motif.name),]
  temp = as.data.frame(cbind(linked.TF$pval.adj,enriched$p.adjust,DATF$DATF.combined.pval.adj))
  temp$V4 = ifelse(sqrt(-log(temp$V1))>20, 20,  sqrt(-log(temp$V1)))
  temp$V5 = ifelse(sqrt(-log(temp$V2))>20, 20,  sqrt(-log(temp$V2))) 
  temp$V6 = ifelse(sqrt(-log(temp$V3))>20, 20,  sqrt(-log(temp$V3)))
  temp=temp[,-c(1:3)]
  names(temp) <- c("Linked TF", "Enriched TF", "DATF")
  
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor=NULL, ...)
  { usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  filter=!is.na(y) & !is.na(x) & !is.nan(x) & !is.nan(y) & !is.infinite(x) & !is.infinite(y)
  r <- cor(x[filter], y[filter])
  r.temp=abs(r)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, 'r = ', txt)
  cex.cor <- 1
  text(0.5, 0.5, txt, cex = 1.2)}
  

pdf("pdf/five/figure4.pdf",width = 5, height = 5)  
  pairs(temp,
        upper.panel = panel.cor, lower.panel = function(x,y){smoothScatter(x,y,add=T)},
        main=('Correlation between TFs'))
  
  dev.off()
  

  linked.TF <- read.csv('linkage_analysis_output_wm/linked.TF.csv')
  linked.gene <- read.csv('linkage_analysis_output_wm/linked.gene.csv')
  linked.peak <- read.table('linkage_analysis_output_wm/linked.peak.output.txt', header = T)
  enriched <- read.csv('linkage_analysis_output_wm/enriched.motif.csv')
  
  head(combined.linked.gene[-1,],20)$gene
    linked.gene.2 <- read.csv('../../HIV2/linkage_analysis_output_wm/linked.gene.csv')
  
  print(xtable(linked.gene[linked.gene$gene %in%   head(combined.linked.gene[-1,],20)$gene, ], digits = c(4,4,4,-4,-4) ,floating=F), include.rownames=FALSE)
  print(xtable(linked.gene.2[linked.gene.2$gene %in%   head(combined.linked.gene[-1,],20)$gene, ], digits = c(4,4,4,-4,-4) ,floating=F), include.rownames=FALSE)
  linked.gene.2=linked.gene.2[match(linked.gene$gene,linked.gene.2$gene),]
  print(xtable(linked.gene.2[linked.gene.2$gene %in%   head(combined.linked.gene[-1,],20)$gene, ], digits = c(4,4,4,-4,-4) ,floating=F), include.rownames=FALSE)
  
  
  linked.TF.2 <- read.csv('../../HIV2/linkage_analysis_output_wm/linked.TF.csv')

  print(xtable(linked.TF[linked.TF$TF %in%   head(combined.linked.TF,20)$TF, ], digits = c(4,4,4,4,-4,-4) ,floating=F), include.rownames=FALSE)
  print(xtable(linked.TF.2[linked.TF.2$TF %in%   head(combined.linked.TF[-1,],20)$TF, ], digits = c(4,4,4,-4,-4) ,floating=F), include.rownames=FALSE)
  linked.TF.2=linked.TF.2[match(linked.TF$TF,linked.TF.2$TF),]
  print(xtable(linked.TF.2[linked.TF.2$TF %in%   head(combined.linked.TF,20)$TF, ], digits = c(4,4,4,4,-4,-4) ,floating=F), include.rownames=FALSE)
  
  
  linked.gene.2 <- read.csv('../../HIV2/linkage_analysis_output_wm/linked.gene.csv')

  
  linked.peak.2 <- read.csv('../../HIV2/linkage_analysis_output_wm/linked.peak.output.txt',header = T)
  linked.peak.2 = linked.peak.2[match(linked.peak$peak,linked.peak.2$peak),]
  
  enriched.2 <- read.csv('../../HIV2/linkage_analysis_output_wm/enriched.motif.csv')
  enriched.2 = enriched.2[match(enriched$motif,enriched.2$motif),]
  
  print(xtable(enriched[enriched$motif.name %in%   head(combined.enriched,20)$TF,c(8,1:7,9) ], digits = c(4,4,4,4,4,4,4,4,-4,-4) ,floating=F), include.rownames=FALSE)
  print(xtable(enriched.2[enriched.2$motif.name %in%   head(combined.enriched,20)$TF,c(8,1:7,9) ], digits = c(4,4,4,4,4,4,4,4,-4,-4) ,floating=F), include.rownames=FALSE)

  
  

  
  p.vector = cbind(linked.gene$pval,linked.gene.2$pval)
  p.vector[p.vector==1]=max(0.99, 1-1/dim(p.vector)[2])
  d <- dim(p.vector)[2]
  linked.gene$cauchy.p.val = NA
  for (i in 1:dim(p.vector)[1]) {
    Sd <- sum(tan((0.5-p.vector[i,])*pi)/d)
    p.global <- pcauchy(Sd,lower.tail = F)
    linked.gene[i,5] <- min(p.global,1)
  }
  
  linked.gene$cauchy.pval.adj = qvalue(linked.gene$cauchy.p.val)$qvalues

  combined.linked.gene = cbind(linked.gene[,c(1,2)],linked.gene.2$cor,linked.gene$cauchy.p.val,linked.gene$cauchy.pval.adj)
  names(combined.linked.gene) <- c("gene", "donor1 cor","donor2 cor", "pval", "pval.adj")
  combined.linked.gene=combined.linked.gene[order(combined.linked.gene$pval),]
  write.csv(combined.linked.gene, file='linkage_analysis_output_wm/combined.linked.gene.csv', row.names = F)
  print(xtable(head(combined.linked.gene[-1,],20) , digits = c(4,4,4,4,-4,-4) ,floating=F), include.rownames=FALSE)
  
  
  gene= head(combined.linked.gene[-1,],20)$gene
  for(i in gene[i]){
    
  }
  

  
  p.vector = cbind(linked.TF$pval,linked.TF.2$pval)
  p.vector[p.vector==1]=max(0.99, 1-1/dim(p.vector)[2])
  d <- dim(p.vector)[2]
  linked.TF$cauchy.p.val = NA
  for (i in 1:dim(p.vector)[1]) {
    Sd <- sum(tan((0.5-p.vector[i,])*pi)/d)
    p.global <- pcauchy(Sd,lower.tail = F)
    linked.TF[i,6] <- min(p.global,1)
  }
  
  linked.TF$cauchy.pval.adj = qvalue(linked.TF$cauchy.p.val)$qvalues
  combined.linked.TF = cbind(linked.TF[,c(1,2,3)],linked.TF.2$cor,linked.TF$cauchy.p.val,linked.TF$cauchy.pval.adj)
  names(combined.linked.TF) <- c("TF", "Motif", "Donor1_cor", "Donor2_cor", "Combined.pval", "Combined.pval.adj")
  combined.linked.TF=combined.linked.TF[order(combined.linked.TF$Combined.pval),]
  write.csv(combined.linked.TF, file='linkage_analysis_output_wm/combined.linked.TF.csv', row.names = F)
  
  print(xtable(head(combined.linked.TF,20) , digits = c(4,4,4,4,4,-4,-4) ,floating=F), include.rownames=FALSE)
  
  
  
  
  p.vector = cbind(linked.TF$pval,linked.TF.2$pval)
  p.vector[p.vector==1]=max(0.99, 1-1/dim(p.vector)[2])
  d <- dim(p.vector)[2]
  linked.TF$cauchy.p.val = NA
  Sd.vec = rep(NA,dim(p.vector)[1])
  for (i in 1:dim(p.vector)[1]) {
    Sd <- sum(tan((0.5-p.vector[i,])*pi)/d)
    Sd.vec[i]=Sd
    p.global <- pcauchy(Sd,lower.tail = F)
    linked.TF[i,6] <- min(p.global,1)
  }
  
  linked.TF$cauchy.pval.adj = qvalue(linked.TF$cauchy.p.val)$qvalues
  combined.linked.TF = cbind(linked.TF[,c(1,2,3)],linked.TF.2$cor,linked.TF$cauchy.p.val,linked.TF$cauchy.pval.adj)
  names(combined.linked.TF) <- c("TF", "Motif", "Donor1_cor", "Donor2_cor", "Combined.pval", "Combined.pval.adj")
  combined.linked.TF=combined.linked.TF[order(combined.linked.TF$Combined.pval),]
  write.csv(combined.linked.TF, file='linkage_analysis_output_wm/combined.linked.TF.csv', row.names = F)
  
  print(xtable(head(combined.linked.TF,20) , digits = c(4,4,4,4,4,-4,-4) ,floating=F), include.rownames=FALSE)
  
  
  
  
  p.vector = cbind(enriched$pval,enriched.2$pval)
  p.vector[p.vector>0.9999]=max(0.9999, 1-1/dim(p.vector)[2])
  d <- dim(p.vector)[2]
  enriched$cauchy.p.val = NA
  for (i in 1:dim(p.vector)[1]) {
    Sd <- sum(tan((0.5-p.vector[i,])*pi)/d)
    p.global <- pcauchy(Sd,lower.tail = F)
    enriched[i,10] <- min(p.global,1)
  }
  
  combined.enriched= as.data.frame(cbind(enriched[,8], enriched$motif, enriched$observed,enriched$background,enriched$percent.observed, enriched$percent.background,
                                         enriched$fold.enrichment, enriched.2$observed,enriched.2$background,enriched.2$percent.observed, enriched.2$percent.background,
                                         enriched.2$fold.enrichment,enriched$cauchy.p.val))
  
  names(combined.enriched) <- c("TF", "Motif", "Donor1_observed", "Donor1_background", "Donor1_%observed", "Donor1_%background","Donor1_enrichment",
                                "Donor2_observed", "Donor2_background", "Donor2_%observed", "Donor2_%background","Donor2_enrichment","Combined.pval")
  
  combined.enriched[,c(3:13)] <- apply(combined.enriched[,c(3:13)] , 2,            # Specify own function within apply
                                      function(x) as.numeric(as.character(x)))
  
  combined.enriched=combined.enriched[order(combined.enriched$Combined.pval),]
  write.csv(combined.enriched, file='linkage_analysis_output_wm/combined.enriched.csv', row.names = F)
  
  print(xtable(head(combined.enriched,20) , digits = c(4,4,0,0,0,4,4,4,0,0,4,4,4,-4) ,floating=F), include.rownames=FALSE)
  
  
  
  
  
 
  p.vector = cbind(DATF$SAHA.pval, DATF$iBET151.pval,DATF$Prostratin.pval,DATF.2$SAHA.pval, DATF.2$iBET151.pval,DATF.2$Prostratin.pval)
  p.vector[p.vector>0.9999]=max(0.9999, 1-1/dim(p.vector)[2])
  d <- dim(p.vector)[2]
  DATF$cauchy.p.val = NA
  for (i in 1:dim(p.vector)[1]) {
    Sd <- sum(tan((0.5-p.vector[i,])*pi)/d)
    p.global <- pcauchy(Sd,lower.tail = F)
    DATF[i,15] <- min(p.global,1)
  }
  DATF$cauchy.p.val.adj = qvalue(DATF$cauchy.p.val,pi0 = 1)$qvalues
    
  DATF.combined =as.data.frame(cbind(DATF$TF,DATF$motif,DATF$iBET151.pval,DATF$SAHA.pval,DATF$Prostratin.pval,DATF.2$iBET151.pval,
                       DATF.2$SAHA.pval,DATF.2$Prostratin.pval,DATF$cauchy.p.val,DATF$cauchy.p.val.adj))

  names(DATF.combined) <- c("TF", "Motif","Donor1_iBET151.pval", "Donor1_SAHA.pval","Donor1_Prostratin.pval",
                            "Donor2_iBET151.pval", "Donor2_SAHA.pval","Donor2_Prostratin.pval", "Combined.pval","Combined.pval.adj")
  DATF.combined[,c(3:10)] <- apply(DATF.combined[,c(3:10)] , 2,            # Specify own function within apply
                                      function(x) as.numeric(as.character(x)))
  DATF.combined=DATF.combined[order(DATF.combined$Combined.pval),]
  write.csv(DATF.combined, file='linkage_analysis_output_wm/DATF.combined.csv', row.names = F)
  
  
  intersect(intersect(combined.linked.TF[combined.linked.TF$Combined.pval.adj<0.01,1],combined.enriched[combined.enriched$Combined.pval<0.05,1]),
            DATF[DATF.combined$Combined.pval.adj<0.01,1])
  
  motif = c("RELA","STAT1", "CREB1", "BACH2", "NFE2" ,"ATF3","JDP2","ATF7","IKZF1", "ZNF75D")
  combined.linked.TF = combined.linked.TF[combined.linked.TF$TF %in% motif,]
  print(xtable(combined.linked.TF[order(combined.linked.TF$Motif),], digits = -4,floating=F), include.rownames=FALSE)
  combined.enriched=combined.enriched[combined.enriched$TF %in% motif,]
  print(xtable(combined.enriched[order(combined.enriched$Motif),], digits = c(4,4,4,4,4,4,4,-4),floating=F), include.rownames=FALSE)
  DATF.combined=DATF.combined[DATF.combined$TF %in% motif,]
  print(xtable(DATF.combined[order(DATF.combined$Motif),] , digits = -4 ,floating=F), include.rownames=FALSE)
  
  
setwd("~/Dropbox/HIV")
setwd("C:/Users/yuchaoj/Dropbox/HIV")

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

load('processed_data/hiv_predict_DE.rda')

# DATF
head(DATF)
colnames(DATF)[3:(ncol(DATF)-2)]=paste0('DATF.', colnames(DATF)[3:(ncol(DATF)-2)])

# Linked TF motif accessibility
linkedTF=read.csv('linkage_analysis_output/linked.TF.csv')
linkedTF=linkedTF[match(DATF$TF,linkedTF$TF),]
colnames(linkedTF)[3:ncol(linkedTF)]=paste0('linkedTF.',colnames(linkedTF)[3:ncol(linkedTF)])

# TF with enriched motif in linked peaks'
rm(enriched.motifs)
enrichedTF=read.csv('linkage_analysis_output/enriched.motif.csv')
enrichedTF=enrichedTF[match(DATF$TF,enrichedTF$motif.name),]
colnames(enrichedTF)[2:(ncol(enrichedTF)-1)]=paste0('enrichedTF.',colnames(enrichedTF)[2:(ncol(enrichedTF)-1)])

TF.output=cbind(DATF, linkedTF, enrichedTF)
TF.output=TF.output[,-c(15,16,20,27)]

alpha=0.1
library(ggvenn)
x=list(linked.TF=linked.TF[linked.TF$pval.adj<=alpha,2],
       enrichedTF=enrichedTF[enrichedTF$enrichedTF.pvalue<=alpha,1],
       DATF=DATF[DATF$DATF.combined.pval.adj<alpha,2])
p1=ggvenn(x, 
          fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)
p1
ggsave(p1, file='output/figure6_venn.pdf', width=3, height=3)

overlapped.motifs=intersect(intersect(linked.TF[linked.TF$pval.adj<=0.1,2], enrichedTF[enrichedTF$enrichedTF.pvalue<=0.1,1]),
                            DATF[DATF$DATF.combined.pval.adj<0.1,2])
overlapped.motifs

linked.TF[match(overlapped.motifs, linked.TF[,2]),]
enrichedTF[match(overlapped.motifs, enrichedTF[,1]),]
DATF[match(overlapped.motifs, DATF[,2]),]

temp=linked.TF[match(overlapped.motifs, linked.TF[,2]),]
temp[,3:5]=signif(temp[,3:5],4)
write.csv(temp, file='output/supp_table4_linkedTF.csv', row.names = F, quote = F)

temp=enrichedTF[match(overlapped.motifs, enrichedTF[,1]),]
temp[,c(4,5,6,7)]=signif(temp[,c(4,5,6,7)],4)
colnames(temp)=gsub('enrichedTF.','', colnames(temp))
colnames(temp)[ncol(temp)]='TF'
temp=temp[,c(ncol(temp),1:(ncol(temp)-1))]
write.csv(temp, file='output/supp_table4_enrichedTF.csv', row.names = F, quote = F)

temp=DATF[match(overlapped.motifs, DATF[,2]),]
temp[,-c(1,2)]=signif(temp[,-c(1,2)],4)
write.csv(temp, file='output/supp_table4_DATF.csv', row.names = F, quote = F)


cauchy.comb=function(Pvals){
  Pvals[Pvals>0.999]=0.999
  is.small<-(Pvals<1e-15)
  Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
  Pvals[is.small]<-1/Pvals[is.small]/pi
  cct.stat<-mean(Pvals)
  pval<-pcauchy(cct.stat,lower.tail = F)
  return(pval)
}

p.cauchy=apply(TF.output[,c("DATF.combined.pval","linkedTF.pval","enrichedTF.pvalue" )],1, cauchy.comb)
p.cauchy.adj=qvalue(p.cauchy)$qvalues

TF.output=cbind(TF.output, p.cauchy, p.cauchy.adj)
TF.output[,-c(1:2)]=signif(TF.output[,-c(1:2)],4)
write.csv(TF.output, file='TF.output.csv', row.names = F)



####################################
#### TF network analysis
####################################

# Release the treshold min_cells, argument passed down to vst in sctransform
hiv <- SCTransform(hiv, verbose = FALSE, return.only.var.genes = FALSE, min_cells=1) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')


# Focus on TFs called by at least two testing
TF.sel=TF.output$TF[order(TF.output$p.cauchy.adj)[1:200]]
length(TF.sel)
TF.sel[grep('SNA', TF.sel)]
TF.sel[grep('CTCF', TF.sel)]

TF.sel=intersect(TF.sel, rownames(hiv@assays$SCT))

# TF exp
TF.exp=hiv@assays$SCT@data[match(TF.sel, rownames(hiv@assays$SCT)),]
# TF motif score
TF.access=hiv@assays$chromvar@data[match(TF.sel, TF),]
length(TF.sel)
dim(TF.exp)
dim(TF.access)

TF.sel=TF.sel[1:50]
TF.exp=TF.exp[1:50,]
TF.access=TF.access[1:50,]

TF.exp.pc=TF.exp%*%svd(TF.exp)$v[,1:10]
TF.access.pc=TF.access%*%svd(TF.access)$v[,1:10]

TF.exp.access.pc=cbind(TF.exp.pc, TF.access.pc)
TF.consensus.pc=as.matrix(TF.exp.access.pc%*%svd(TF.exp.access.pc)$v)
TF.cor=cor(t(TF.consensus.pc))

library(pheatmap)
pdf(file='output/figure6_pheatmap.pdf', width=8, height=8)
pheatmap(TF.cor)
dev.off()

# igraph plotting: https://kateto.net/wp-content/uploads/2016/06/Polnet%202016%20R%20Network%20Visualization%20Workshop.pdf
# adjacency matrix based on correlations: https://davetang.org/muse/2017/03/16/matrix-to-adjacency-list-in-r/

# replace upper triangle of the matrix
# with number that can be used for filtering
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
set.seed(2)
plot(net, layout = layout_components(net), edge.width = E(net)$weight)
pdf(file='output/figure6_network.pdf', width=10, height=10)
plot(net, layout = layout_components(net), edge.width = E(net)$weight, vertex.shape="none")
dev.off()


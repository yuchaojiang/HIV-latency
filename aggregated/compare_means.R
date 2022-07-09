library(tidyverse)
library(ggplot2)
library(patchwork)
library(Seurat)
library(SeuratObject)
library(Signac)
# install.packages("rstatix")
# library(rstatix)
# install.packages("ggpubr")
# library(ggpubr)

file.seurat <- "derived_data/hiv_predict.rds"
# file.vln <- "figures/vln_hiv_cond_anova.pdf"
file.anova <- "derived_data/anova.csv"
file.tukey <- "derived_data/tukey.csv"

seurat <- file.seurat %>% readRDS()
d <- seurat@meta.data %>%
    select(pRNA.sqrt, group)

# get one-way ANOVA p-value
# anova <- aov(pRNA.sqrt ~ group, data = d)
lm <- lm(pRNA.sqrt ~ group, data = d)
anova <- anova(lm)
file.anova %>% write.csv(anova, file = .)

# get Tukey multiple comparisons of means
tukey <- aov(pRNA.sqrt ~ group, data = d) %>%
    TukeyHSD()
file.tukey %>% write.csv(tukey$group, file = .)

# https://www.biostars.org/p/458261/
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
# https://github.com/kassambara/ggpubr/issues/102

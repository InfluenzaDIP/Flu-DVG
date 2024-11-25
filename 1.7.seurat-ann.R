library(Seurat)
library(patchwork)
library(ggpubr)
library(tidyverse)
library(magrittr)
library(future)

set.seed(1515)
options(future.globals.maxSize = 80000 * 1024^2, seed = TRUE)
plan("multicore", workers = 40)
plan()

theme.2 = theme(plot.margin=unit(c(0,0,0,0), 'cm'),
                plot.title = element_blank(),
                text=element_text(size = 8),
                axis.line = element_blank(),
                axis.text.y = element_blank(), 
                axis.ticks.y = element_blank(), 
                axis.text.x = element_blank(), 
                axis.ticks.x = element_blank(),
                panel.border = element_rect(fill=NA, color="black", size=0.2, linetype="solid"),
                legend.key.size = unit(0.3,"cm"),
                legend.margin = unit(0.3,"cm"))
theme.4 = theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), 'cm'),
                plot.title = element_blank(),
                text=element_text(size = 8),
                axis.text.y = element_text(size = 8), 
                # axis.ticks.y = element_blank(), 
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
                # axis.ticks.x = element_blank(),
                axis.title = element_blank(),
                panel.grid = element_blank(),
                legend.key.size = unit(0.3,"cm"),
                legend.margin = unit(0.3,"cm"),
                legend.title=element_blank())


# 1. Load data
combined <- readRDS("./results/1.5.1.merged-singleR-20231102.rds")
ann = read.csv("./results/1.7.1.ann.tsv", sep = '\t')
dim(ann)

col.pal.ann = read.csv('./results/1.7.1.color.tsv', sep = '\t', header = F)
col.use.ann = col.pal.ann[, 2]
names(col.use.ann) = col.pal.ann[, 1]



# 2.ann
combined@meta.data$Ann <- ann[match(combined@meta.data$RNA_snn_res.1.2, ann$Cluster), 2]
combined@meta.data$Ann = factor(combined@meta.data$Ann, levels = col.pal.ann[, 1])
head(combined@meta.data)
table(combined@meta.data$Ann)
saveRDS(combined, file = "./results/1.7.2.merged-ann-20231102.rds")
write.csv(combined@meta.data, file = "./results/1.7.2.meta.csv")



# 3.plot
pdf("./results/01-07-01-UMAP-harmony-ann.pdf", 
    width = 20/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "Ann", 
    # label = TRUE, 
    pt.size = 0.2, label.size = 2.5, repel = FALSE) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    scale_color_manual(values = col.use.ann) + 
    theme.2
dev.off()

data = as.data.frame(table(combined@meta.data$Ann))
colnames(data) = c('ann', 'Count')
pdf("./results/01-07-02-bar.ann.withlegend.pdf", 
    width = 20/2.54, height = 8/2.54)
ggplot(data, aes(ann, Count, fill = ann))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = col.use.ann) + 
    theme_bw()+
    theme.4
dev.off()

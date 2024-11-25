library(Seurat)
library(patchwork)
library(tidyverse)
library(magrittr)
library(future)

set.seed(1515)
options(future.globals.maxSize = 80000 * 1024^2, seed = TRUE)
plan("multicore", workers = 40)
plan()


# 1. Load data
combined <- readRDS("./results/1.1.2.merged-raw-20231102.rds")
dim(combined)


# 2. plot
theme.1 = theme(plot.margin=unit(c(0,0,0,0), 'cm'),
                plot.title = element_text(size = 8),
                text = element_text(size = 8),
                axis.text = element_text(size = 6), 
                axis.title = element_blank(),
                panel.grid = element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")
pdf("./results/01-02-01-violin-raw.pdf", 
    width = 20/2.54, height = 5/2.54)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(combined, group.by = "sample", pt.size = 0,
                       features = plot.featrures[i]) + theme.1 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1)     
violin
dev.off()


# 3. qc
minGene=500
maxGene=5000
minUMI=1000
maxUMI=50000
pctMT=10
### Êý¾ÝÖÊ¿Ø
combined <- subset(combined, 
                   subset = nCount_RNA > minUMI & 
                   nCount_RNA < maxUMI & 
                   nFeature_RNA > minGene & 
                   # nFeature_RNA < maxGene & 
                   percent.mt < pctMT)
dim(combined)


# 4. plot
pdf("./results/01-02-02-violin-qc.pdf", 
    width = 20/2.54, height = 5/2.54)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(combined, group.by = "sample", pt.size = 0,
                       features = plot.featrures[i]) + theme.1 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1)     
violin
dev.off()


# 5. save
saveRDS(combined, file = "./results/1.2.1.merged-qc-20231102.rds")

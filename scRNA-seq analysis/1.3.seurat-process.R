library(Seurat)
library(homologene)
library(patchwork)
library(ggpubr)
library(tidyverse)
library(magrittr)
library(future)

set.seed(1515)
options(future.globals.maxSize = 80000 * 1024^2, seed = TRUE)
plan("multicore", workers = 40)
plan()


# 1. Load data
combined <- readRDS("/workdir/wangph/proj/20231102_Jida/SC-3sample/1.2.1.merged-qc-20231102.rds")
dim(combined)

genes = read.csv('/workdir/wangph/proj/20230320_TB_resister/4.seurat/v4-norm-harmony/1.overall/1.11.1.genes.unwanted.txt', sep = '\t', header = F)
genes = genes$V1
# genes = genes[genes$V2 != 'Dissociation', 1]
genes = homologene(genes, inTax = 9606, outTax = 10090)$'10090'
combined = combined[!rownames(combined) %in% genes, ]
dim(combined)


# 2. 
combined <- NormalizeData(combined) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
# combined <- SC-3sampleTransform(combined, conserve.memory=TRUE)


# 3. pca
combined <- RunPCA(combined, verbose = F)
saveRDS(combined, file = "/workdir/wangph/proj/20231102_Jida/SC-3sample/1.3.1.merged-pca-20231102.rds")

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-01-elbow.pdf", 
    width = 10/2.54, height = 10/2.54)
ElbowPlot(combined, ndims = 50)
dev.off()


# 4.
pc.num=1:30
combined <- RunTSNE(combined, dims=pc.num) %>% 
    RunUMAP(dims=pc.num) %>% 
    FindNeighbors(dims = pc.num) %>%
    FindClusters(resolution = c(0.1, 0.2, 0.5, 0.8, 1, 1.2, 1.5, 2))
saveRDS(combined, file = "/workdir/wangph/proj/20231102_Jida/SC-3sample/1.3.2.merged-umap-20231102.rds")


# 5. plot
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

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-02-UMAP-cluster0.1.pdf", 
    width = 12/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "RNA_snn_res.0.2", 
    label = TRUE, pt.size = 0.2, label.size = 2.5, repel = FALSE) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    theme.2
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-02-UMAP-cluster0.2.pdf", 
    width = 12/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "RNA_snn_res.0.2", 
    label = TRUE, pt.size = 0.2, label.size = 2.5, repel = FALSE) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    theme.2
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-02-UMAP-cluster0.5.pdf", 
    width = 12/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "RNA_snn_res.0.5", 
    label = TRUE, pt.size = 0.2, label.size = 2.5, repel = FALSE) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    theme.2
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-02-UMAP-cluster0.8.pdf", 
    width = 12/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "RNA_snn_res.0.8", 
    label = TRUE, pt.size = 0.2, label.size = 2.5, repel = FALSE) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    theme.2
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-02-UMAP-cluster1.0.pdf", 
    width = 12/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "RNA_snn_res.1", 
    label = TRUE, pt.size = 0.2, label.size = 2.5, repel = FALSE) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    theme.2
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-02-UMAP-cluster1.2.pdf", 
    width = 12/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "RNA_snn_res.1.2", 
    label = TRUE, pt.size = 0.2, label.size = 2.5, repel = FALSE) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    theme.2
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-02-UMAP-cluster1.5.pdf", 
    width = 12/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "RNA_snn_res.1.5", 
    label = TRUE, pt.size = 0.2, label.size = 2.5, repel = FALSE) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    theme.2
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-02-UMAP-cluster2.0.pdf", 
    width = 12/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "RNA_snn_res.2", 
    label = TRUE, pt.size = 0.2, label.size = 2.5, repel = FALSE) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    theme.2
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-03-UMAP-group.pdf", 
    width = 32/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "RNA_snn_res.1", split.by = "sample", 
    label = FALSE, ncol = 4) + 
    theme.2
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-04-feature-cellcycle.pdf", 
    width = 16/2.54, height = 8/2.54)
# feature = c('PCNA', 'MKI67')
feature = c('Pcna', 'Mki67')
plots = list()
for(i in seq_along(feature)){
    plots[[i]] = FeaturePlot(combined, features = feature[i], label = FALSE, 
                             pt.size = 0.2, label.size = 2, repel = FALSE) + 
    NoLegend() +
    theme.2 +
    theme(axis.title = element_blank(), 
          plot.title = element_text(size = 8))
}
wrap_plots(plots = plots, nrow=1)
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-05-feature-celltype.pdf", 
    width = 24/2.54, height = 24/2.54)
feature = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
feature = homologene(feature, inTax = 9606, outTax = 10090)$'10090'
plots = list()
for(i in seq_along(feature)){
    plots[[i]] = FeaturePlot(combined, features = feature[i], label = FALSE, 
                             pt.size = 0.2, label.size = 2, repel = FALSE) + 
    NoLegend() +
    theme.2 +
    theme(axis.title = element_blank(), 
          plot.title = element_text(size = 8))
}
wrap_plots(plots = plots, nrow=3)
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-05-feature-celltype-b.pdf", 
    width = 16/2.54, height = 16/2.54)
feature = c('MS4A1', 'CD27', 'CD38', 'IGHG1')
feature = homologene(feature, inTax = 9606, outTax = 10090)$'10090'
plots = list()
for(i in seq_along(feature)){
    plots[[i]] = FeaturePlot(combined, features = feature[i], label = FALSE, 
                             pt.size = 0.2, label.size = 2, repel = FALSE) + 
    NoLegend() +
    theme.2 +
    theme(axis.title = element_blank(), 
          plot.title = element_text(size = 8))
}
wrap_plots(plots = plots, nrow=2)
dev.off()


# 6. markers
feature = read.table('/workdir/wangph/proj/20231102_Jida/SC-3sample/markers-epi.txt')
feature = feature$V1
pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-06-marker-epi.pdf", 
    width = 40/2.54, height = 40/2.54)
plots = list()
for(i in seq_along(feature)){
    plots[[i]] = FeaturePlot(combined, features = feature[i], label = FALSE, 
                             pt.size = 0.2, label.size = 2, repel = FALSE) + 
    NoLegend() +
    theme.2 +
    theme(axis.title = element_blank(), 
          plot.title = element_text(size = 8))
}
wrap_plots(plots = plots, nrow=5)
dev.off()

feature = read.table('/workdir/wangph/proj/20231102_Jida/SC-3sample/markers-fib.txt')
feature = feature$V1
pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-06-marker-fib.pdf", 
    width = 24/2.54, height = 24/2.54)
plots = list()
for(i in seq_along(feature)){
    plots[[i]] = FeaturePlot(combined, features = feature[i], label = FALSE, 
                             pt.size = 0.2, label.size = 2, repel = FALSE) + 
    NoLegend() +
    theme.2 +
    theme(axis.title = element_blank(), 
          plot.title = element_text(size = 8))
}
wrap_plots(plots = plots, nrow=3)
dev.off()

feature = read.table('/workdir/wangph/proj/20231102_Jida/SC-3sample/markers-imm.txt')
feature = feature$V1
pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-06-marker-imm.pdf", 
    width = 48/2.54, height = 40/2.54)
plots = list()
for(i in seq_along(feature)){
    plots[[i]] = FeaturePlot(combined, features = feature[i], label = FALSE, 
                             pt.size = 0.2, label.size = 2, repel = FALSE) + 
    NoLegend() +
    theme.2 +
    theme(axis.title = element_blank(), 
          plot.title = element_text(size = 8))
}
wrap_plots(plots = plots, nrow=5)
dev.off()

feature = read.table('/workdir/wangph/proj/20231102_Jida/SC-3sample/markers-neu.txt')
feature = feature$V1
pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-03-06-marker-neu.pdf", 
    width = 24/2.54, height = 24/2.54)
plots = list()
for(i in seq_along(feature)){
    plots[[i]] = FeaturePlot(combined, features = feature[i], label = FALSE, 
                             pt.size = 0.2, label.size = 2, repel = FALSE) + 
    NoLegend() +
    theme.2 +
    theme(axis.title = element_blank(), 
          plot.title = element_text(size = 8))
}
wrap_plots(plots = plots, nrow=3)
dev.off()

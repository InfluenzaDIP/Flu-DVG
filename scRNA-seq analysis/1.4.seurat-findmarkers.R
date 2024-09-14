library(Seurat)
library(harmony)
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
combined <- readRDS("/workdir/wangph/proj/20231102_Jida/SC-3sample/1.3.2.merged-umap-20231102.rds")
dim(combined)


# 2. 
Idents(combined) <- combined@meta.data$RNA_snn_res.0.5
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.markers, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-04-01-markers-res0.5-20231102.csv')

Idents(combined) <- combined@meta.data$RNA_snn_res.1
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.markers, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-04-01-markers-res1.0-20231102.csv')

Idents(combined) <- combined@meta.data$RNA_snn_res.2
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.markers, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-04-01-markers-res2.0-20231102.csv')


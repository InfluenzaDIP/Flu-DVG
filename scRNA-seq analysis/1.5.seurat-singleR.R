library(Seurat)
library(SingleR)
library(celldex)
library(pheatmap)
library(dplyr)
library(ggpubr)
library(future)

set.seed(1515)
options(future.globals.maxSize = 80000 * 1024^2, seed = TRUE)
plan("multicore", workers = 40)
plan()


# 1. Load data
combined <- readRDS("/workdir/wangph/proj/20231102_Jida/SC-3sample/1.3.2.merged-umap-20231102.rds")
dim(combined)


# 2. 
be.se <- BlueprintEncodeData()
dice.se <- DatabaseImmuneCellExpressionData()
hpca.se <- HumanPrimaryCellAtlasData()
mi.se <- MonacoImmuneData()
nh.se <- NovershternHematopoieticData()
mouseImmu <- ImmGenData()
mouseRNA <- MouseRNAseqData()


# 3. 进行singleR注释
aggr_for_SingleR <- GetAssayData(combined, slot = "data") ##获取标准化矩阵
aggr.singler <- SingleR(test = aggr_for_SingleR, 
                        ref = list(mouseImmu = mouseImmu, mouseRNA = mouseRNA), 
                        labels = list(mouseImmu$label.main, mouseRNA$label.main)) 
combined@meta.data$singleR.main = aggr.singler$labels
aggr.singler <- SingleR(test = aggr_for_SingleR, 
                        ref = list(mouseImmu = mouseImmu, mouseRNA = mouseRNA), 
                        labels = list(mouseImmu$label.fine, mouseRNA$label.fine)) 
combined@meta.data$singleR.fine = aggr.singler$labels
head(combined@meta.data)
# saveRDS(combined, file = "/workdir/wangph/proj/20231102_Jida/SC-3sample/1.5.1.merged-singleR-20231102.rds")


# 4. 整理注释
combined@meta.data$singleR.main[combined@meta.data$singleR.main == 'B cells, pro'] = 'B cells'
combined@meta.data$singleR.main[combined@meta.data$singleR.main == 'DC'] = 'Dendritic cells'
combined@meta.data$singleR.main[combined@meta.data$singleR.main == 'ILC'] = 'NK cells'
combined@meta.data$singleR.main[combined@meta.data$singleR.main == 'NKT'] = 'T cells'
combined@meta.data$singleR.main[combined@meta.data$singleR.main == 'Tgd'] = 'T cells'

combined@meta.data$singleR.fine = gsub(" \\(.*\\)","",combined@meta.data$singleR.fine)
combined@meta.data$singleR.fine = gsub(" activated","",combined@meta.data$singleR.fine)
combined@meta.data$singleR.fine[combined@meta.data$singleR.fine == 'aNSCs'] = 'Neural progenitor cells'
combined@meta.data$singleR.fine[combined@meta.data$singleR.fine == 'DC'] = 'Dendritic cells'
combined@meta.data$singleR.fine[combined@meta.data$singleR.fine == 'Ependymal'] = 'Ependymal cells'
combined@meta.data$singleR.fine[combined@meta.data$singleR.fine == 'Fibroblasts senescent'] = 'Fibroblasts'
combined@meta.data$singleR.fine[combined@meta.data$singleR.fine == 'ILC'] = 'NK cells'
combined@meta.data$singleR.fine[combined@meta.data$singleR.fine == 'NKT'] = 'T cells'
combined@meta.data$singleR.fine[combined@meta.data$singleR.fine == 'NPCs'] = 'Neural progenitor cells'
combined@meta.data$singleR.fine[combined@meta.data$singleR.fine == 'OPCs'] = 'Oligodendrocyte progenitor cells'
combined@meta.data$singleR.fine[combined@meta.data$singleR.fine == 'qNSCs'] = 'Neural progenitor cells'
combined@meta.data$singleR.fine[combined@meta.data$singleR.fine == 'Tgd'] = 'T cells'

saveRDS(combined, file = "/workdir/wangph/proj/20231102_Jida/SC-3sample/1.5.1.merged-singleR-20231102.rds")


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

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-05-01-UMAP-singleR.main.pdf", 
    width = 13.5/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "singleR.main", 
    label = TRUE, pt.size = 0.2, label.size = 2.5, repel = TRUE) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    theme.2
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-05-01-UMAP-singleR.fine.pdf", 
    width = 15/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "singleR.fine", 
    label = TRUE, pt.size = 0.2, label.size = 2.5, repel = TRUE) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    theme.2
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-05-02-UMAP-singleR.main-group.pdf", 
    width = 28/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "singleR.main", split.by = "sample", 
    label = FALSE, ncol = 3) + 
    theme.2
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-05-02-UMAP-singleR.fine-group.pdf", 
    width = 28/2.54, height = 8/2.54)
DimPlot(combined, reduction = "umap", 
    group.by = "singleR.fine", split.by = "sample", 
    label = FALSE, ncol = 3) + 
    theme.2
dev.off()

# 6. save table
tbl = table(combined@meta.data$singleR.main, combined@meta.data$RNA_snn_res.2)
write.csv(tbl, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-05-03-singleR-harmony-cluster2.0.main.csv', 
          quote = TRUE)
tbl = table(combined@meta.data$singleR.main, combined@meta.data$RNA_snn_res.1)
write.csv(tbl, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-05-03-singleR-harmony-cluster1.0.main.csv', 
          quote = TRUE)
tbl = table(combined@meta.data$singleR.main, combined@meta.data$RNA_snn_res.0.5)
write.csv(tbl, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-05-03-singleR-harmony-cluster0.5.main.csv', 
          quote = TRUE)
tbl = table(combined@meta.data$singleR.fine, combined@meta.data$RNA_snn_res.2)
write.csv(tbl, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-05-03-singleR-harmony-cluster2.0.fine.csv', 
          quote = TRUE)
tbl = table(combined@meta.data$singleR.fine, combined@meta.data$RNA_snn_res.1)
write.csv(tbl, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-05-03-singleR-harmony-cluster1.0.fine.csv', 
          quote = TRUE)
tbl = table(combined@meta.data$singleR.fine, combined@meta.data$RNA_snn_res.0.5)
write.csv(tbl, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-05-03-singleR-harmony-cluster0.5.fine.csv', 
          quote = TRUE)

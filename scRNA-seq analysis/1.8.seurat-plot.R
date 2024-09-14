library(Seurat)
library(ggplot2)
library(ggpointdensity)
library(pheatmap)
library(RColorBrewer)
library(ggsci)
library(paletteer)
library(patchwork)
library(ggpubr)
library(ggalluvial)
library(cowplot)
library(dplyr)
library(ggrepel)
library(msigdbr)
library(GSVA)

set.seed(1515)
options(future.globals.maxSize = 80000 * 1024^2, seed = TRUE)
plan("multicore", workers = 40)
plan()

theme.1 = theme(plot.margin=unit(c(0,0,0,0), 'cm'),
                plot.title = element_text(size = 8),
                text = element_text(size = 8),
                axis.text = element_text(size = 6), 
                axis.title = element_blank(),
                panel.grid = element_blank())
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
theme.3 = theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), 'cm'),
                plot.title = element_blank(),
                text=element_text(size = 8),
                axis.text.y = element_blank(), 
                axis.ticks.y = element_blank(), 
                axis.text.x = element_blank(), 
                axis.ticks.x = element_blank(),
                axis.title = element_blank(),
                panel.grid = element_blank(),
                legend.key.size = unit(0.3,"cm"),
                legend.margin = unit(0.3,"cm"),
                legend.title=element_blank())
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

col_lite = function(cols, ds=0.4, dv=0.7) {
  cols = rgb2hsv(col2rgb(cols))
  cols["v", ] = cols["v", ] + dv*(1 - cols["v", ])
  cols["s", ] = ds*cols["s", ]
  apply(cols, 2, function(x) hsv(x[1], x[2], x[3]))
}

save_pheatmap_pdf <- function(x, filename, width=8, height=12) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


# 1. Load data
combined <- readRDS("/workdir/wangph/proj/20231102_Jida/SC-3sample/1.7.2.merged-ann-20231102.rds")
dim(combined)
combined@meta.data$Virus = factor(combined@meta.data$Virus, levels = c('Yes', 'No'))
head(combined@meta.data)

col.pal.ann = read.csv('/workdir/wangph/proj/20231102_Jida/SC-3sample/1.7.1.color.tsv', sep = '\t', header = F)
col.use.ann = col.pal.ann[, 2]
names(col.use.ann) = col.pal.ann[, 1]


# 2. prop
prop.ann <- as.data.frame(table(as.character(combined@meta.data$Ann), combined@meta.data$sample))
colnames(prop.ann)<-c("cell_type","sample","count")
prop.ann <- transform(prop.ann, percentage = ave(count, sample, FUN = prop.table)) 
prop.ann$cell_type = factor(prop.ann$cell_type, 
                           levels = col.pal.ann[, 1])
head(prop.ann)


prop.virus <- as.data.frame(table(as.character(combined@meta.data$Virus), combined@meta.data$sample))
colnames(prop.virus)<-c("cell_type","sample","count")
prop.virus <- transform(prop.virus, percentage = ave(count, sample, FUN = prop.table)) 
prop.virus$cell_type = factor(prop.virus$cell_type, 
                           levels = c('Yes', 'No'))
head(prop.virus)

count.sample = as.data.frame(table(combined@meta.data$sample))
colnames(count.sample)<-c("sample","count")


# 3. dge
Idents(combined) = combined@meta.data$Ann
markers <- FindAllMarkers(combined, min.pct = 0.5, logfc.threshold = 0)
write.csv(markers, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-00-dge-all-20231102.csv')


# 4. GSVA
hallmark.class = read.csv('/workdir/wangph/proj/20220418_TB/SC/5.seurat/v7-final-v2/1.data/04-hallmark-class-wj01.csv', row.names = 1)
hallmark.class$class = factor(hallmark.class$class, levels = c("immune", "metabolic/hypoxia/ROS", "cell growth and proliferation",
                                                                      "cell development", "DNA damage",  "cellular component"))
hallmark = msigdbr(species = "Mus musculus", category = "H") 
hallmark.df = subset(hallmark, select = c("gs_name","gene_symbol")) %>% as.data.frame()
hallmark.list = split(hallmark.df$gene_symbol, hallmark.df$gs_name)
gobp = msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") 
gobp.df = subset(gobp, select = c("gs_name","gene_symbol")) %>% as.data.frame()
gobp.list = split(gobp.df$gene_symbol, gobp.df$gs_name)
kegg = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG") 
kegg.df = subset(kegg, select = c("gs_name","gene_symbol")) %>% as.data.frame()
kegg.list = split(kegg.df$gene_symbol, kegg.df$gs_name)

expr = AverageExpression(combined, assays = "RNA", slot = "data")[[1]]
expr = expr[rowSums(expr)>0,]  #选取非零基因
expr = as.matrix(expr)
head(expr)

gsva.hallmark = gsva(expr, hallmark.list, method="gsva", kcdf="Gaussian") # "Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
gsva.kegg = gsva(expr, kegg.list, method="gsva", kcdf="Gaussian") 
gsva.gobp = gsva(expr, gobp.list, method="gsva", kcdf="Gaussian")
write.csv(gsva.hallmark, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-00-gsva-hallmark.csv')
write.csv(gsva.kegg, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-00-gsva-kegg.csv')
write.csv(gsva.gobp, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-00-gsva-gobp.csv')


Idents(combined) = combined@meta.data$sample
expr = AverageExpression(combined, assays = "RNA", slot = "data")[[1]]
expr = expr[rowSums(expr)>0,] 
write.csv(expr, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-00-exp1.csv')

combined@meta.data$temp = paste0(combined@meta.data$sample, '_', combined@meta.data$Ann)
Idents(combined) = combined@meta.data$temp
expr = AverageExpression(combined, assays = "RNA", slot = "data")[[1]]
expr = expr[rowSums(expr)>0,] 
write.csv(expr, 
          '/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-00-exp2.csv')



# 5. plot
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")
pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-01-violin-qc.pdf", 
    width = 12/2.54, height = 12/2.54)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(combined, group.by = "sample", split.by =  "sample",  pt.size = 0,
                       features = plot.featrures[i]) +  
                       scale_fill_manual(values = pal_nejm(alpha = 0.7)(9))  + 
                       theme.1  # +
                       # NoLegend()
}
violin <- wrap_plots(plots = plots, nrow=2)     
violin
dev.off()

## UMAP图
pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-02-UMAP.cluster.withlegend.pdf", 
    width = 11/2.54, height = 8/2.54)
DimPlot(combined, group.by = "RNA_snn_res.1.2",reduction = "umap", 
        label = TRUE, 
        pt.size = 0.2, label.size = 2.5, repel = FALSE) + 
    # NoLegend() +
    labs(x = "UMAP1", y = "UMAP2") + 
    # scale_color_manual(values = col.use.cluster) +
    theme.2
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-03-UMAP.ann.withlegend.pdf", 
    width = 19/2.54, height = 8/2.54)
DimPlot(combined, group.by = "Ann", reduction = "umap", 
        label = FALSE, repel = TRUE, 
        pt.size = 0.2, label.size = 2.5) + 
    # NoLegend() +
    labs(x = "UMAP1", y = "UMAP2") + 
    scale_color_manual(values = col.use.ann) +
    theme.2
dev.off()
data = as.data.frame(table(combined@meta.data$Ann))
colnames(data) = c('Ann', 'Count')
pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-03-bar.ann.withlegend.pdf", 
    width = 20/2.54, height = 8/2.54)
ggplot(data, aes(Ann, Count, fill = Ann))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = col.use.ann) +
    theme_bw()+
    theme.4
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-04-UMAP.virus.withlegend.pdf", 
    width = 10/2.54, height = 8/2.54)
DimPlot(combined, group.by = "Virus", reduction = "umap", 
        pt.size = 0.2, label.size = 2.5) + 
    # NoLegend() +
    labs(x = "UMAP1", y = "UMAP2") + 
    scale_color_manual(values = pal_nejm(alpha = 0.7)(8)) +
    theme.2
dev.off()
data = subset(x = combined, subset = Virus == 'Yes')
data = as.data.frame(table(data@meta.data$sample, data@meta.data$Ann))
colnames(data) = c('Sample', 'Ann', 'Count')
pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-04-bar.virus.withlegend.pdf", 
    width = 30/2.54, height = 8/2.54)
ggplot(data, aes(Sample, Count, fill = Ann))+
    geom_bar(stat = "identity", position = 'dodge')+
    geom_text(aes(label = Count), vjust = - 0.2, position = position_dodge(.9), size = 2) +
    scale_fill_manual(values = col.use.ann) +
    # scale_y_log10() + 
    theme_bw()+
    theme.4
dev.off()

## UMAP（分样本）
pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-05-UMAP.ann.seperated_by_group.pdf", 
    width = 22/2.54, height = 8/2.54)
DimPlot(combined, group.by = "Ann", split.by = "sample", reduction = "umap", label = FALSE)  + 
    NoLegend() + 
    scale_color_manual(values = col.use.ann) +
    theme.2
dev.off()
pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-05-UMAP.virus.seperated_by_group.pdf", 
    width = 22/2.54, height = 8/2.54)
DimPlot(combined, group.by = "Virus", split.by = "sample", reduction = "umap", label = FALSE)  + 
    NoLegend() + 
    scale_color_manual(values = pal_nejm(alpha = 0.7)(8)) +
    theme.2
dev.off()

## UMAP密度图（分样本）
data <- cbind(Embeddings(object=combined[['umap']]),FetchData(combined,'sample'))
pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-06-UMAP.density.seperated_by_group.pdf", 
    width = 25/2.54, height = 8/2.54)
ggplot(data = data, mapping = aes(x = UMAP_1,y = UMAP_2)) + 
    geom_pointdensity(size = 0.001) + 
    scale_color_paletteer_c("viridis::plasma") +
    facet_grid(~sample) + 
    theme_minimal() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme.2
dev.off()


# prop
pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-07-flow.prop_by_ann.withlegend.pdf", 
    width = 14/2.54, height = 8/2.54)
ggplot(prop.ann, aes(x = sample, y = percentage, fill = cell_type, stratum=cell_type, alluvium=cell_type)) +
    geom_col(width = 0.5, color="black", size = 0.1) +
    geom_flow(width = 0.5, alpha=0.2, knot.pos=0)+
    scale_fill_manual(values = col.use.ann) +
    scale_x_discrete(labels = paste(count.sample$sample, '\n(', count.sample$count, ')'))+ 
    theme_bw()+
    theme.4
dev.off()

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-07-flow.prop_by_virus.withlegend.pdf", 
    width = 7/2.54, height = 8/2.54)
ggplot(prop.virus, aes(x = sample, y = percentage, fill = cell_type, stratum=cell_type, alluvium=cell_type)) +
    geom_col(width = 0.5, color="black", size = 0.1) +
    geom_flow(width = 0.5, alpha=0.2, knot.pos=0)+
    scale_fill_manual(values = pal_nejm(alpha = 0.7)(8)) +
    scale_x_discrete(labels = paste(count.sample$sample, '\n(', count.sample$count, ')'))+ 
    theme_bw()+
    theme.4
dev.off()

## 气泡图
markers <- read.csv('/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-00-dge-all-20231102.csv', )
top <- markers %>% 
    group_by(cluster) %>% 
    top_n(2, avg_log2FC)
top=top[!duplicated(top$gene),]
top$cluster = factor(top$cluster, levels = col.pal.ann[, 1])
top = top[order(top$cluster),]
# head(top)
select_genes_all=split(top$gene, top$cluster)
# select_genes_all
n <- c()
for(x in select_genes_all){
    n<-c(n, x)
}
pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-08-dot.markers.pdf", 
    width = 24/2.54, height = 10/2.54)
DotPlot(object = combined, group.by = "Ann", features = n, assay = "RNA", dot.scale = 3) + 
    scale_color_distiller(palette = "RdYlBu") + 
    theme.1 +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
dev.off()

## 火山图
markers$cluster = factor(markers$cluster, levels = col.pal.ann[, 1])
markers$label = ifelse(markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) >= 0.5,
                ifelse(markers$avg_log2FC > 0,'red','blue'), 'gray')
color = c(red = "#CB181D", gray = "#BFC0C2", blue = "#2171B5")
marker.list = split(markers, markers$cluster)
for(i in c(1:length(marker.list))){
    sig_dge.pos = marker.list[[i]][marker.list[[i]]$avg_log2FC >= 0.5 & 
                                   marker.list[[i]]$p_val_adj <= 0.05, ] %>% 
                                   slice_max(n = 10, order_by = avg_log2FC)
    sig_dge.neg = marker.list[[i]][marker.list[[i]]$avg_log2FC <= -0.5 & 
                                   marker.list[[i]]$p_val_adj <= 0.05, ] %>% 
                                   slice_min(n = 10, order_by = avg_log2FC)
    sig_dge = rbind(sig_dge.pos, sig_dge.neg)
    # print(sig_dge)
    # print(head(marker.list[[i]]))
    
    pdf(paste("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-09-volcano.deg-", 
                            sub('/', '_', marker.list[[i]][1, 6]),  ".pdf", sep = ""), 
        width = 8/2.54, height = 8/2.54)
    plot = ggplot(marker.list[[i]], aes(avg_log2FC, -log10(p_val_adj), col = label)) +
        geom_point(size = 1) + 
        # xlim(-5, 5) +
        scale_color_manual(values = color) +
        geom_hline(yintercept = -log10(0.05), lty = 4, col = "grey", lwd = 0.6) +
        geom_vline(xintercept = c(-0.5, 0.5), lty = 4, col = "grey", lwd = 0.6) + 
        labs(x = "Log2 (fold change)",y = "-log10 (adjusted p-value)") +
        geom_label_repel(data = sig_dge, aes(label=gene), size = 3, max.overlaps = 30, col = "#3C3C3C") +
        theme_bw()+
        theme(legend.position = "none",
              panel.grid = element_blank(),
              axis.title =  element_text(size = 8),
              axis.text = element_text(size = 8),
              strip.text = element_text(size = 8),
              strip.background = element_blank())
    print(plot)
    dev.off()
}

#labels for top genes
top10.pos = markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
top10.neg = markers %>%
    group_by(cluster) %>%
    slice_min(n = 10, order_by = avg_log2FC)
top10 = rbind(top10.pos, top10.neg)
top10 = top10[top10$p_val_adj <= 0.05,]

#添加X轴的cluster色块标签：
dfcol = data.frame(x = c(1: length(col.pal.ann[, 1])), 
	y = 0, 
	label = col.pal.ann[, 1])
mycol = col.use.ann

pdf("/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-10-volcano.deg.pdf", 
    width = 100/2.54, height = 16/2.54)
ggplot()+
    geom_jitter(data = markers, aes(x = cluster, y = avg_log2FC, color = label),
                size = 0.85, width = 0.4)+
    geom_text_repel(data = top10, aes(x = cluster, y = avg_log2FC, label = gene),
                    force = 1.2, max.overlaps = 100,
                    arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last"))+
    scale_color_manual(name = NULL, values = c(red = "#CB181D", gray = "#BFC0C2", blue = "#2171B5"))+ 
    geom_tile(data = dfcol, aes(x = x,y = y), 
              height = 0.4, color = "black", fill = mycol, alpha = 1, show.legend = F)+
    labs(x = "Cluster",y = "average logFC")+
    geom_text(data = dfcol,
              aes(x = x, y = y, label = label),
              size = 4, color ="Black") +
    theme_minimal()+
    theme(axis.title = element_text(size = 12, color = "black", face = "bold"),
          axis.line.y = element_line(color = "black", size = 1),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none",
          legend.direction = "vertical",
          legend.justification = c(1,0),
          legend.text = element_text(size = 15)
     )
dev.off()

## GSVA
breaksList = seq(-0.8, 0.8, by = 0.01)
colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
p = pheatmap(gsva.hallmark[rownames(hallmark.class), ], 
             scale = "none",
             na_col = "white", 
             color = colours, 
             breaks = breaksList,
             cluster_cols = FALSE, 
             cluster_row = FALSE, 
             annotation_row = hallmark.class,
             # labels_col = paste(count.celltype$celltype, '(', count.celltype$count, ')'),
             # display_numbers = as.matrix(ifelse(data_wide_adjp <= 0.05, "*", "")), 
             border=FALSE)
save_pheatmap_pdf(p, "/workdir/wangph/proj/20231102_Jida/SC-3sample/01-08-11-gsva.pdf",
	12, 12)

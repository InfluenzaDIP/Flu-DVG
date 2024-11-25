library(Seurat)
library(msigdbr)
library(GSVA)
library(patchwork)
library(ggpubr)
library(tidyverse)
library(magrittr)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(future)


set.seed(1515)
options(future.globals.maxSize = 80000 * 1024^2, seed = TRUE)
plan("multicore", workers = 40)
plan()

save_pheatmap_pdf <- function(x, filename, width=8, height=12) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


# 1. Load data
combined = readRDS("./results/1.7.2.merged-ann-20231102.rds")
dim(combined)
comparisons = read.csv('./results/1.9.0.comparison-2.txt', sep = '\t', header = FALSE)
head(comparisons)
col.pal.ann = read.csv('./results/1.7.1.color.tsv', sep = '\t', header = F)
col.use.ann = col.pal.ann[, 2]
names(col.use.ann) = col.pal.ann[, 1]


# 2. findmarkers
Idents(combined) = combined@meta.data$Ann
for(i in c(1:nrow(comparisons))){
	for (cluster in col.pal.ann[, 1]){
		sub = subset(x = combined, subset = Ann == cluster)
		markers = FindMarkers(sub, 
			ident.1 = comparisons$V1[i], 
			ident.2 = comparisons$V2[i],
			group.by = 'sample', 
			min.pct = 0, 
			logfc.threshold = 0)
		write.csv(markers, 
			paste("./results/01-09-1", i, "-dge-", comparisons$V1[i], 'vs', comparisons$V2[i], '-', cluster, ".csv", sep = ""))
	}
}


# 3. plot
color = c(red = "#CB181D", gray = "#BFC0C2", blue = "#2171B5")
deg_all = c()
data_all = c()
for(i in c(1:nrow(comparisons))){
	deg.pos = c()
	deg.neg = c()
	data = c()
	for (cluster in col.pal.ann[, 1]){
		markers = read.csv(paste("./results/01-09-1", i, "-dge-", comparisons$V1[i], 'vs', comparisons$V2[i], '-', cluster, ".csv", sep = ""),
	                                   header = T, row.names = 1)
		markers$gene = rownames(markers)
		markers$cluster = cluster
		markers$label = ifelse(markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) >= 0.5,
                                ifelse(markers$avg_log2FC > 0,'red','blue'), 'gray')
		data = rbind(data, markers)
		sig_dge.pos = markers[markers$avg_log2FC >= 0.5 & 
                                      markers$p_val_adj <= 0.05, ] %>% 
                                      slice_max(n = 10, order_by = avg_log2FC)
		sig_dge.neg = markers[markers$avg_log2FC <= -0.5 & 
                                      markers$p_val_adj <= 0.05, ] %>% 
                                      slice_min(n = 10, order_by = avg_log2FC)
		sig_dge = rbind(sig_dge.pos, sig_dge.neg)
		deg.pos = c(deg.pos, sig_dge.pos$gene)
		deg.neg = c(deg.neg, sig_dge.neg$gene)

		pdf(paste("./results/01-09-1", i, "-dge-", comparisons$V1[i], 'vs', comparisons$V2[i], '-', cluster, ".pdf", sep = ""), 
			width = 10/2.54, height = 10/2.54)
		plot = ggplot(markers, aes(avg_log2FC, -log10(p_val_adj), col = label)) +
                       geom_point(size = 1) + 
                       # xlim(-5, 5) +
                       scale_color_manual(values = color) +
                       geom_hline(yintercept = -log10(0.05), lty = 4, col = "grey", lwd = 0.6) +
                       geom_vline(xintercept = c(-0.5, 0.5), lty = 4, col = "grey", lwd = 0.6) + 
                       labs(x = "Log2 (fold change)", y = "-log10 (adjusted p-value)") +
                       geom_label_repel(data = sig_dge, aes(label = gene), size = 3, max.overlaps = 30, col = "#3C3C3C") +
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


	features = unique(deg.pos)
	data_long_log2fc = data[data$gene %in% features, c('cluster', 'gene', 'avg_log2FC')]
	data_long_log2fc$cluster =  factor(data_long_log2fc$cluster, levels = col.pal.ann[, 1])
	data_wide_log2fc = spread(data_long_log2fc, cluster, avg_log2FC)
	# data_wide_log2fc
	rownames(data_wide_log2fc) = data_wide_log2fc$gene
	data_wide_log2fc = data_wide_log2fc[features, c(2:ncol(data_wide_log2fc))]
	data_long_adjp = data[data$gene %in% features, c('cluster', 'gene', 'p_val_adj')]
	data_long_adjp$cluster =  factor(data_long_adjp$cluster, levels = col.pal.ann[, 1])
	data_wide_adjp = spread(data_long_adjp, cluster, p_val_adj)
	# data_wide_adjp
	rownames(data_wide_adjp) = data_wide_adjp$gene
	data_wide_adjp = data_wide_adjp[features, c(2:ncol(data_wide_adjp))]

	breaksList = seq(-3, 3, by = 0.01)
	colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
	p1 = pheatmap(data_wide_log2fc, 
             	scale = "none",
             	na_col = "white", 
             	color = colours, 
             	breaks = breaksList,
             	cluster_cols = FALSE, 
             	# cluster_row = FALSE, 
             	# labels_col = paste(count.celltype$celltype, '(', count.celltype$count, ')'),
             	display_numbers = as.matrix(ifelse(data_wide_adjp <= 0.05, "*", "")), 
             	border=FALSE)
	save_pheatmap_pdf(p1, paste("./results/01-09-2", i, "-dge-", comparisons$V1[i], 'vs', comparisons$V2[i], ".pos.pdf", sep = ""), 
		20/2.54, 30/2.54)


	features = unique(deg.neg)
	data_long_log2fc = data[data$gene %in% features, c('cluster', 'gene', 'avg_log2FC')]
	data_long_log2fc$cluster =  factor(data_long_log2fc$cluster, levels = col.pal.ann[, 1])
	data_wide_log2fc = spread(data_long_log2fc, cluster, avg_log2FC)
	# data_wide_log2fc
	rownames(data_wide_log2fc) = data_wide_log2fc$gene
	data_wide_log2fc = data_wide_log2fc[features, c(2:ncol(data_wide_log2fc))]
	data_long_adjp = data[data$gene %in% features, c('cluster', 'gene', 'p_val_adj')]
	data_long_adjp$cluster =  factor(data_long_adjp$cluster, levels = col.pal.ann[, 1])
	data_wide_adjp = spread(data_long_adjp, cluster, p_val_adj)
	# data_wide_adjp
	rownames(data_wide_adjp) = data_wide_adjp$gene
	data_wide_adjp = data_wide_adjp[features, c(2:ncol(data_wide_adjp))]

	breaksList = seq(-3, 3, by = 0.01)
	colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
	p2 = pheatmap(data_wide_log2fc, 
             	scale = "none",
             	na_col = "white", 
             	color = colours, 
             	breaks = breaksList,
             	cluster_cols = FALSE, 
             	# cluster_row = FALSE, 
             	# labels_col = paste(count.celltype$celltype, '(', count.celltype$count, ')'),
             	display_numbers = as.matrix(ifelse(data_wide_adjp <= 0.05, "*", "")), 
             	border=FALSE)
	save_pheatmap_pdf(p2, paste("./results/01-09-2", i, "-dge-", comparisons$V1[i], 'vs', comparisons$V2[i], ".neg.pdf", sep = ""), 
		20/2.54, 30/2.54)


	features = unique(c(deg.pos, deg.neg))
	data_long_log2fc = data[data$gene %in% features, c('cluster', 'gene', 'avg_log2FC')]
	data_long_log2fc$cluster =  factor(data_long_log2fc$cluster, levels = col.pal.ann[, 1])
	data_wide_log2fc = spread(data_long_log2fc, cluster, avg_log2FC)
	# data_wide_log2fc
	rownames(data_wide_log2fc) = data_wide_log2fc$gene
	data_wide_log2fc = data_wide_log2fc[features, c(2:ncol(data_wide_log2fc))]
	data_long_adjp = data[data$gene %in% features, c('cluster', 'gene', 'p_val_adj')]
	data_long_adjp$cluster =  factor(data_long_adjp$cluster, levels = col.pal.ann[, 1])
	data_wide_adjp = spread(data_long_adjp, cluster, p_val_adj)
	# data_wide_adjp
	rownames(data_wide_adjp) = data_wide_adjp$gene
	data_wide_adjp = data_wide_adjp[features, c(2:ncol(data_wide_adjp))]

	breaksList = seq(-3, 3, by = 0.01)
	colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
	p3 = pheatmap(data_wide_log2fc, 
             	scale = "none",
             	na_col = "white", 
             	color = colours, 
             	breaks = breaksList,
             	cluster_cols = FALSE, 
             	# cluster_row = FALSE, 
             	# labels_col = paste(count.celltype$celltype, '(', count.celltype$count, ')'),
             	display_numbers = as.matrix(ifelse(data_wide_adjp <= 0.05, "*", "")), 
             	border=FALSE)
	save_pheatmap_pdf(p3, paste("./results/01-09-2", i, "-dge-", comparisons$V1[i], 'vs', comparisons$V2[i], ".pdf", sep = ""), 
		20/2.54, 50/2.54)

	data$comparison = paste0(comparisons$V1[i], 'vs', comparisons$V2[i])
	data_all = rbind(data_all, data)
	deg_all = c(deg_all, features)
	write.csv(data, paste("./results/01-09-2", i, "-dge-", comparisons$V1[i], 'vs', comparisons$V2[i], ".csv", sep = ""))
}

data_all$comparison = paste0(data_all$cluster, ' (', data_all$comparison, ')')
write.csv(data_all, './results/01-09-24-deg-all.csv')
data_all = read.csv('./results/01-09-24-deg-all.csv', row.names = 1)


data_long_log2fc = data_all[data_all$gene %in% features, c('comparison', 'gene', 'avg_log2FC')]
data_long_log2fc$comparison =  factor(data_long_log2fc$comparison, levels = unique(data_all$comparison))
data_wide_log2fc = spread(data_long_log2fc, comparison, avg_log2FC)
# data_wide_log2fc
rownames(data_wide_log2fc) = data_wide_log2fc$gene
data_wide_log2fc = data_wide_log2fc[features, c(2:ncol(data_wide_log2fc))]
data_long_adjp = data_all[data_all$gene %in% features, c('comparison', 'gene', 'p_val_adj')]
data_long_adjp$comparison =  factor(data_long_adjp$comparison, levels = unique(data_all$comparison))
data_wide_adjp = spread(data_long_adjp, comparison, p_val_adj)
# data_wide_adjp
rownames(data_wide_adjp) = data_wide_adjp$gene
data_wide_adjp = data_wide_adjp[features, c(2:ncol(data_wide_adjp))]

breaksList = seq(-3, 3, by = 0.01)
colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
p = pheatmap(data_wide_log2fc, 
	scale = "none",
	na_col = "white", 
	color = colours, 
	breaks = breaksList,
	cluster_cols = FALSE, 
	# cluster_row = FALSE, 
	# labels_col = paste(count.celltype$celltype, '(', count.celltype$count, ')'),
	display_numbers = as.matrix(ifelse(data_wide_adjp <= 0.05, "*", "")), 
	border=FALSE,
	gaps_col = c(26, 52))
save_pheatmap_pdf(p, "./results/01-09-24-dge-BvsC.pdf", 
	50/2.54, 60/2.54)

features = unique(deg_all)
# data_all$comparison = paste0(data_all$cluster, ' (', data_all$comparison, ')')
data_long_log2fc = data_all[data_all$gene %in% features, c('comparison', 'gene', 'avg_log2FC')]
data_long_log2fc$comparison =  factor(data_long_log2fc$comparison, levels = unique(data_all$comparison))
data_wide_log2fc = spread(data_long_log2fc, comparison, avg_log2FC)
# data_wide_log2fc
rownames(data_wide_log2fc) = data_wide_log2fc$gene
data_wide_log2fc = data_wide_log2fc[features, c(2:ncol(data_wide_log2fc))]
data_long_adjp = data_all[data_all$gene %in% features, c('comparison', 'gene', 'p_val_adj')]
data_long_adjp$comparison =  factor(data_long_adjp$comparison, levels = unique(data_all$comparison))
data_wide_adjp = spread(data_long_adjp, comparison, p_val_adj)
# data_wide_adjp
rownames(data_wide_adjp) = data_wide_adjp$gene
data_wide_adjp = data_wide_adjp[features, c(2:ncol(data_wide_adjp))]

breaksList = seq(-3, 3, by = 0.01)
colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
p = pheatmap(data_wide_log2fc, 
	scale = "none",
	na_col = "white", 
	color = colours, 
	breaks = breaksList,
	cluster_cols = FALSE, 
	# cluster_row = FALSE, 
	# labels_col = paste(count.celltype$celltype, '(', count.celltype$count, ')'),
	display_numbers = as.matrix(ifelse(data_wide_adjp <= 0.05, "*", "")), 
	border=FALSE,
	gaps_col = c(26, 52))
save_pheatmap_pdf(p, "./results/01-09-24-dge-all.pdf", 
	50/2.54, 100/2.54)

features =  read.table('./results/01-09-10-GeneSelected.txt')
features = features$V1
# data_all$comparison = paste0(data_all$cluster, ' (', data_all$comparison, ')')
data_long_log2fc = data_all[data_all$gene %in% features, c('comparison', 'gene', 'avg_log2FC')]
data_long_log2fc$comparison =  factor(data_long_log2fc$comparison, levels = unique(data_all$comparison))
data_wide_log2fc = spread(data_long_log2fc, comparison, avg_log2FC)
# data_wide_log2fc
rownames(data_wide_log2fc) = data_wide_log2fc$gene
data_wide_log2fc = data_wide_log2fc[features, c(2:ncol(data_wide_log2fc))]
data_long_adjp = data_all[data_all$gene %in% features, c('comparison', 'gene', 'p_val_adj')]
data_long_adjp$comparison =  factor(data_long_adjp$comparison, levels = unique(data_all$comparison))
data_wide_adjp = spread(data_long_adjp, comparison, p_val_adj)
# data_wide_adjp
rownames(data_wide_adjp) = data_wide_adjp$gene
data_wide_adjp = data_wide_adjp[features, c(2:ncol(data_wide_adjp))]

breaksList = seq(-2, 2, by = 0.01)
colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
p = pheatmap(data_wide_log2fc, 
	scale = "none",
	na_col = "white", 
	color = colours, 
	breaks = breaksList,
	cluster_cols = FALSE, 
	cluster_row = FALSE, 
	# labels_col = paste(count.celltype$celltype, '(', count.celltype$count, ')'),
	display_numbers = as.matrix(ifelse(data_wide_adjp <= 0.05, "*", "")), 
	border=FALSE,
	gaps_col = c(26, 52))
save_pheatmap_pdf(p, "./results/01-09-24-dge-selecetd.pdf", 
	50/2.54, 40/2.54)

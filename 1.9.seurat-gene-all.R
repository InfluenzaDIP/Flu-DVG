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


# 2. findmarkers
# Idents(combined) = combined@meta.data$ann
for(i in c(1:nrow(comparisons))){	
	markers = FindMarkers(combined, 
		ident.1 = comparisons$V1[i], 
		ident.2 = comparisons$V2[i],
		group.by = 'sample', 
		min.pct = 0, 
		logfc.threshold = 0)
	write.csv(markers, 
		paste("./results/01-09-01-dge-", comparisons$V1[i], 'vs', comparisons$V2[i], ".csv", sep = ""))
}


# 3. plot
color = c(red = "#CB181D", gray = "#BFC0C2", blue = "#2171B5")
for(i in c(1:nrow(comparisons))){
	markers = read.csv(paste("./results/01-09-01-dge-", comparisons$V1[i], 'vs', comparisons$V2[i], ".csv", sep = ""),
                                   header = T, row.names = 1)
	markers$gene = rownames(markers)
	markers$label = ifelse(markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) >= 0.5,
                                ifelse(markers$avg_log2FC > 0,'red','blue'), 'gray')
	sig_dge.pos = markers[markers$avg_log2FC >= 0.5 & 
                                      markers$p_val_adj <= 0.05, ] %>% 
                                      slice_max(n = 10, order_by = avg_log2FC)
	sig_dge.neg = markers[markers$avg_log2FC <= -0.5 & 
                                      markers$p_val_adj <= 0.05, ] %>% 
                                      slice_min(n = 10, order_by = avg_log2FC)
	sig_dge = rbind(sig_dge.pos, sig_dge.neg)
	pdf(paste("./results/01-09-01-dge-", comparisons$V1[i], 'vs', comparisons$V2[i], ".pdf", sep = ""), 
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

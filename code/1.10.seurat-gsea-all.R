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
library(clusterProfiler)
library(openxlsx)
library(cowplot)


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
comparisons = read.csv('./results/1.9.0.comparison-2.txt', sep = '\t', header = FALSE)
head(comparisons)
col.pal.ann = read.csv('./results/1.7.1.color.tsv', sep = '\t', header = F)
col.use.ann = col.pal.ann[, 2]
names(col.use.ann) = col.pal.ann[, 1]

hallmark.class = read.csv('/workdir/wangph/proj/20220418_TB/SC/5.seurat/v7-final-v2/1.data/04-hallmark-class-wj01.csv', row.names = 1)
hallmark.class$class = factor(hallmark.class$class, levels = c("immune", "metabolic/hypoxia/ROS", "cell growth and proliferation",
                                                                      "cell development", "DNA damage",  "cellular component"))
hallmark <- msigdbr(species = "Mus musculus", category = "H") 
hallmark.df <- subset(hallmark, select = c("gs_name","gene_symbol")) %>% as.data.frame()
hallmark.list <- split(hallmark.df$gene_symbol, hallmark.df$gs_name)

gobp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") 
gobp.df <- subset(gobp, select = c("gs_name","gene_symbol")) %>% as.data.frame()
gobp.list <- split(gobp.df$gene_symbol, gobp.df$gs_name)

kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG") 
kegg.df <- subset(kegg, select = c("gs_name","gene_symbol")) %>% as.data.frame()
kegg.list <- split(kegg.df$gene_symbol, kegg.df$gs_name)


# 2. GSEA
for(i in c(1:nrow(comparisons))){
	data = read.csv(paste0("./results/01-09-0", i, "-dge-", comparisons$V1[i], 'vs', comparisons$V2[i], ".csv"), row.names = 1)
	data = data[data$avg_log2FC != 0, ]
	data = data[order(data$avg_log2FC, decreasing = T),]
	genelist <- structure(data$avg_log2FC, names = rownames(data))
	hallmark.gsea <- GSEA(genelist, 
	                      TERM2GENE = hallmark.df, 
	                      eps = 0,
	                      pvalueCutoff  = 1)
	write.csv(data.frame(hallmark.gsea), 
	          file = paste0("./results/02-10-4", i, 
	                         "-gsea-hallmark-", comparisons$V1[i], 'vs', comparisons$V2[i], ".csv"),
	                         row.names = FALSE)

	gobp.gsea <- GSEA(genelist, 
	                  TERM2GENE = gobp.df, 
	                  eps = 0,
	                  pvalueCutoff  = 1)
	write.csv(data.frame(gobp.gsea), 
	          file = paste0("./results/02-10-4", i, 
	                         "-gsea-gobp-", comparisons$V1[i], 'vs', comparisons$V2[i], ".csv"),
	                         row.names = FALSE)

	kegg.gsea <- GSEA(genelist, 
                          TERM2GENE = kegg.df, 
                          eps = 0,
                          pvalueCutoff  = 1)
	write.csv(data.frame(kegg.gsea), 
	          file = paste0("./results/02-10-4", i, 
	                         "-gsea-kegg-", comparisons$V1[i], 'vs', comparisons$V2[i], ".csv"),
	                         row.names = FALSE)
}


# 3.Hallmark plot
hallmark.BvsA = read.csv("./results/02-10-41-gsea-hallmark-BvsA.csv")
hallmark.CvsA = read.csv("./results/02-10-42-gsea-hallmark-CvsA.csv")
hallmark.BvsC = read.csv("./results/02-10-43-gsea-hallmark-BvsC.csv")

# create a matric of NES 
hallmark.nes.dflist = merge(hallmark.BvsA[, c('ID', 'NES')], hallmark.CvsA[, c('ID', 'NES')], 
                            by='ID', all.x =TRUE, all.y=TRUE)
hallmark.nes.dflist = merge(hallmark.nes.dflist, hallmark.BvsC[, c('ID', 'NES')], 
                            by='ID', all.x =TRUE, all.y=TRUE)
rownames(hallmark.nes.dflist) = hallmark.nes.dflist$ID
colnames(hallmark.nes.dflist)[2:ncol(hallmark.nes.dflist)] = c('BvsA', 'CvsA', 'BvsC')
hallmark.nes.dflist = hallmark.nes.dflist[, 2:ncol(hallmark.nes.dflist)]
hallmark.nes.dflist[is.na(hallmark.nes.dflist)] = 0

# create a matric of p.adjust 
hallmark.adjp.dflist = merge(hallmark.BvsA[, c('ID', 'p.adjust')], hallmark.CvsA[, c('ID', 'p.adjust')], 
                            by='ID', all.x =TRUE, all.y=TRUE)
hallmark.adjp.dflist = merge(hallmark.adjp.dflist, hallmark.BvsC[, c('ID', 'p.adjust')], 
                            by='ID', all.x =TRUE, all.y=TRUE)
rownames(hallmark.adjp.dflist) = hallmark.adjp.dflist$ID
colnames(hallmark.adjp.dflist)[2:ncol(hallmark.adjp.dflist)] = c('BvsA', 'CvsA', 'BvsC')
hallmark.adjp.dflist = hallmark.adjp.dflist[, 2:ncol(hallmark.adjp.dflist)]
hallmark.adjp.dflist[is.na(hallmark.adjp.dflist)] = 1


hallmark.class = read.csv('/workdir/wangph/proj/20220418_TB/SC/5.seurat/v7-final-v2/1.data/04-hallmark-class-wj01.csv')
hallmark.class = hallmark.class[hallmark.class$Hallmark %in% rownames(hallmark.nes.dflist), ]
n = hallmark.class[, 1]
hallmark.class = as.data.frame(hallmark.class[, 2])
rownames(hallmark.class) = n
colnames(hallmark.class) = 'class'
hallmark.class$class = factor(hallmark.class$class, levels = c("immune", "metabolic/hypoxia/ROS", "cell growth and proliferation",
                                                               "cell development", "DNA damage",  "cellular component"))

th = max(max(hallmark.nes.dflist), abs(min(hallmark.nes.dflist))) + 0.5
breaksList = seq(-th, th, by = 0.0001)
colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
p = pheatmap(hallmark.nes.dflist[rownames(hallmark.class), ], 
             scale = "none",
             na_col = "white", 
             color = colours, 
             breaks = breaksList,
             cluster_cols = FALSE, 
             cluster_row = FALSE, 
             annotation_row = hallmark.class, 
             # labels_col = paste(count.celltype$celltype, '(', count.celltype$count, ')'),
             display_numbers = as.matrix(ifelse(hallmark.adjp.dflist[rownames(hallmark.class), ] <= 0.05, "*", "")), 
             border=FALSE)
save_pheatmap_pdf(p, 
                  "./results/02-10-44-gsea-hallmark.pdf", 
                  8, 12)





















# 3.Hallmark plot
for(i in c(1:nrow(comparisons))){
	# hallmark
	path <- paste0("./results/02-10-0", i, 
                       "-gsea-hallmark-", comparisons$V1[i], 'vs', comparisons$V2[i], ".xlsx")
	# getting data from sheets
	sheets <- getSheetNames(path)
	hallmark.gsea.dflist <- lapply(sheets, read.xlsx, xlsxFile=path)
	# assigning names to data frame
	names(hallmark.gsea.dflist) <- col.pal.ann[, 1]
	# get all significant gene sets
	hallmark.sig.list <- lapply(hallmark.gsea.dflist, FUN = function(x){
    	                            res = x[x$p.adjust <= 0.05, 'ID']
    	                            return(res)
	})
	# merge all significant gene sets and remove the duplicated
	hallmark.uniq.list = c()
	for(j in seq_along(hallmark.sig.list)){
	    hallmark.uniq.list <- c(hallmark.uniq.list, hallmark.sig.list[[j]])
	}
	hallmark.uniq.list = sort(unique(hallmark.uniq.list))
	# create a matric of NES 
	hallmark.nes.dflist = hallmark.gsea.dflist[[1]][, c('ID', 'NES')]
	for(j in c(2:length(hallmark.gsea.dflist))){
	    hallmark.nes.dflist <- merge(hallmark.nes.dflist, hallmark.gsea.dflist[[j]][, c('ID', 'NES')], by='ID', all.x =TRUE, all.y=TRUE)
	}
	rownames(hallmark.nes.dflist) = hallmark.nes.dflist$ID
	colnames(hallmark.nes.dflist)[2:ncol(hallmark.nes.dflist)] = col.pal.ann[, 1]
	hallmark.nes.dflist = hallmark.nes.dflist[, 2:ncol(hallmark.nes.dflist)]
	# hallmark.nes.dflist = hallmark.nes.dflist[hallmark.uniq.list, 2:ncol(hallmark.nes.dflist)]
	# create a matric of p.adjust 
	hallmark.adjp.dflist = hallmark.gsea.dflist[[1]][, c('ID', 'p.adjust')]
	for(j in c(2:length(hallmark.gsea.dflist))){
	    hallmark.adjp.dflist <- merge(hallmark.adjp.dflist, hallmark.gsea.dflist[[j]][, c('ID', 'p.adjust')], by='ID', all.x =TRUE, all.y=TRUE)
	}
	rownames(hallmark.adjp.dflist) = hallmark.adjp.dflist$ID
	colnames(hallmark.adjp.dflist)[2:ncol(hallmark.adjp.dflist)] = col.pal.ann[, 1]
	hallmark.adjp.dflist = hallmark.adjp.dflist[, 2:ncol(hallmark.adjp.dflist)]
	# hallmark.adjp.dflist = hallmark.adjp.dflist[hallmark.uniq.list, 2:ncol(hallmark.adjp.dflist)]
	hallmark.adjp.dflist[is.na(hallmark.adjp.dflist)] = 1

	save(list = c("hallmark.nes.dflist", "hallmark.adjp.dflist"), 
             file = paste0("./results/02-10-1", i, 
                       "-gsea-hallmark-", comparisons$V1[i], 'vs', comparisons$V2[i], ".Rdata"))
}

for(i in c(1:nrow(comparisons))){
	# hallmark
	load(paste0("./results/02-10-1", i, 
	            "-gsea-hallmark-", comparisons$V1[i], 'vs', comparisons$V2[i], ".Rdata"))

	hallmark.class = read.csv('/workdir/wangph/proj/20220418_TB/SC/5.seurat/v7-final-v2/1.data/04-hallmark-class-wj01.csv')
	hallmark.class = hallmark.class[hallmark.class$Hallmark %in% rownames(hallmark.nes.dflist), ]
	n = hallmark.class[, 1]
	hallmark.class = as.data.frame(hallmark.class[, 2])
	rownames(hallmark.class) = n
	colnames(hallmark.class) = 'class'
	hallmark.class$class = factor(hallmark.class$class, levels = c("immune", "metabolic/hypoxia/ROS", "cell growth and proliferation",
                                                               "cell development", "DNA damage",  "cellular component"))

	th = max(max(hallmark.nes.dflist), abs(min(hallmark.nes.dflist))) + 0.5
	breaksList = seq(-th, th, by = 0.0001)
	colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
	p = pheatmap(hallmark.nes.dflist[rownames(hallmark.class), ], 
             	scale = "none",
             	na_col = "white", 
             	color = colours, 
             	breaks = breaksList,
             	cluster_cols = FALSE, 
             	cluster_row = FALSE, 
             	annotation_row = hallmark.class, 
             	# labels_col = paste(count.celltype$celltype, '(', count.celltype$count, ')'),
             	display_numbers = as.matrix(ifelse(hallmark.adjp.dflist[rownames(hallmark.class), ] <= 0.05, "*", "")), 
             	border=FALSE)
	save_pheatmap_pdf(p, 
                          paste0("./results/02-10-2", i, 
	                  "-gsea-hallmark-", comparisons$V1[i], 'vs', comparisons$V2[i], ".pdf"), 
                          12, 12)

}

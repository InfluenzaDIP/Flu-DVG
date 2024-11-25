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
	data = read.csv(paste0("./results/01-09-2", i, "-dge-", comparisons$V1[i], 'vs', comparisons$V2[i], ".csv"), row.names = 1)
	data = data[data$avg_log2FC != 0, ]
	data$cluster = factor(data$cluster, levels = col.pal.ann[, 1])
	deg.list <- split(data, f = data$cluster)

	# hallmark
	hallmark.gsea.list <- lapply(deg.list, FUN = function(x){
                                     x = x[order(x$avg_log2FC, decreasing = T),]
                             	     genelist <- structure(x$avg_log2FC, names = x$gene)
                             	     res <- GSEA(genelist, 
                                                 TERM2GENE = hallmark.df, 
                		                 eps = 0,
                		                 pvalueCutoff  = 1)
                                     return(res) #按差异倍数进行排序并返回GSEA结果
                               })
	wb <- createWorkbook()
	for(j in seq_along(hallmark.gsea.list)){
		res <- data.frame(hallmark.gsea.list[[j]])
		if (nchar(names(hallmark.gsea.list)[j]) <= 13){
			sheetName = names(hallmark.gsea.list)[j]
		}else{
			sheetName = substr(names(hallmark.gsea.list)[j], 1, 13)
		}
		addWorksheet(wb, sheetName = sheetName)
		writeData(wb, sheetName, res, rowNames = FALSE)
	}
	saveWorkbook(wb, 
                     file = paste0("./results/02-10-0", i, 
                                   "-gsea-hallmark-", comparisons$V1[i], 'vs', comparisons$V2[i], ".xlsx"), 
                     overwrite = TRUE)

	# GOBP
	gobp.gsea.list <- lapply(deg.list, FUN = function(x){
                         	x = x[order(x$avg_log2FC, decreasing = T),]
                         	genelist <- structure(x$avg_log2FC, names = x$gene)
                         	res <- GSEA(genelist, 
                                     	TERM2GENE = gobp.df, 
                                     	eps = 0,
                                     	pvalueCutoff  = 1)
                         	return(res) #按差异倍数进行排序并返回GSEA结果
                   	})
	wb <- createWorkbook()
	for(j in seq_along(gobp.gsea.list)){
		res <- data.frame(gobp.gsea.list[[j]])
		if (nchar(names(gobp.gsea.list)[j]) <= 13){
			sheetName = names(gobp.gsea.list)[j]
		}else{
			sheetName = substr(names(gobp.gsea.list)[j], 1, 13)
		}
		addWorksheet(wb, sheetName = sheetName)
		writeData(wb, sheetName, res, rowNames = FALSE)
		}
	saveWorkbook(wb, 
                     file = paste0("./results/02-10-0", i, 
                                   "-gsea-gobp-", comparisons$V1[i], 'vs', comparisons$V2[i], ".xlsx"), 
                     overwrite = TRUE)

	# KEGG
	kegg.gsea.list <- lapply(deg.list, FUN = function(x){
                         	x = x[order(x$avg_log2FC, decreasing = T),]
                         	genelist <- structure(x$avg_log2FC, names = x$gene)
                         	res <- GSEA(genelist, 
                                     	TERM2GENE = kegg.df, 
                                     	eps = 0,
                                     	pvalueCutoff  = 1)
                         	return(res) #按差异倍数进行排序并返回GSEA结果
                         	})
	wb <- createWorkbook()
	for(j in seq_along(kegg.gsea.list)){
		res <- data.frame(kegg.gsea.list[[j]])
		if (nchar(names(kegg.gsea.list)[j]) <= 13){
			sheetName = names(kegg.gsea.list)[j]
		}else{
			sheetName = substr(names(kegg.gsea.list)[j], 1, 13)
		}
		addWorksheet(wb, sheetName = sheetName)
		writeData(wb, sheetName, res, rowNames = FALSE)
		}
	saveWorkbook(wb, 
                     file = paste0("./results/02-10-0", i, 
                                   "-gsea-kegg-", comparisons$V1[i], 'vs', comparisons$V2[i], ".xlsx"), 
                     overwrite = TRUE)
}



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


# 4. GO plot
sets = read.table('./results/02-10-00-gobp.txt', sep = '\t')

data.gobp.all = c()
for(i in c(1:nrow(comparisons))){
	# gobp
	path <- paste0("./results/02-10-0", i, 
                       "-gsea-gobp-", comparisons$V1[i], 'vs', comparisons$V2[i], ".xlsx")
	# getting data from sheets
	sheets <- getSheetNames(path)
	gobp.gsea.dflist <- lapply(sheets, read.xlsx, xlsxFile=path)
	# assigning names to data frame
	names(gobp.gsea.dflist) <- col.pal.ann[, 1]
	# for(j in seq_along(gobp.gsea.dflist)){
	for(j in c(2, 3, 4, 7, 10, 11)){
		print(names(gobp.gsea.dflist)[j])
		selected = sets[sets$V2 == names(gobp.gsea.dflist)[j], 1]
		data = gobp.gsea.dflist[[j]]
		data = data[data$ID %in% selected, c('ID', 'NES', 'p.adjust')]
		data$comparison = paste0(comparisons$V1[i], 'vs', comparisons$V2[i])
		data$CellType = names(gobp.gsea.dflist)[j]
		# print(head(data, 2))
		data.gobp.all = rbind(data.gobp.all, data)
	}
}

for(j in c(2, 3, 4, 7, 10, 11)){
	data = data.gobp.all[data.gobp.all$CellType == names(gobp.gsea.dflist)[j], ]
	# print(head(data, 2))
	print(names(gobp.gsea.dflist)[j])

	data$ID = factor(data$ID, levels = rev(sets[sets$V2 == names(gobp.gsea.dflist)[j], 1]))
	data$comparison = factor(data$comparison, levels = c('BvsA', 'CvsA', 'BvsC'))
	g1 = ggplot(data, aes(comparison, ID)) +
	     geom_point(aes(size = -log10(p.adjust), color = p.adjust <= 0.05, fill = NES), shape = 21, stroke = 1) +
	     scale_y_discrete(label = substring(levels(data$ID), 0, 60)) +
             scale_color_manual(values = c('#FFFFFF00', 'black')) +
             scale_fill_gradient2(# limits = c(-3.2, 3.2),
                                  midpoint = 0, 
                                  low = 'navy', 
                                  mid = 'white',
                                  high = 'firebrick3') + 

	     theme_cowplot()
	pdf(paste0("./results/02-10-31-gsea-gobp-dotplot-", 
	           names(gobp.gsea.dflist)[j], ".pdf"), 
	    width = 12, height = 10)
	print(g1)
	dev.off()

	g2 = ggplot(data, aes(-log10(p.adjust), ID, shape = comparison)) +
	     geom_point(aes(color = NES, size = 2)) +
	     geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'lightgray') + 
	     scale_y_discrete(label = substring(levels(data$ID), 0, 60)) +
             scale_color_gradient2(# limits = c(-3.2, 3.2),
                                  midpoint = 0, 
                                  low = 'navy', 
                                  mid = 'white',
                                  high = 'firebrick3') +
	     theme_cowplot()
	pdf(paste0("./results/02-10-32-gsea-gobp-dotplot-", 
	           names(gobp.gsea.dflist)[j], ".pdf"), 
	    width = 12, height = 10)
	print(g2)
	dev.off()

	data$ID = factor(data$ID, levels = rev(sets[sets$V2 == names(gobp.gsea.dflist)[j], 1]))
	data$comparison = factor(data$comparison, levels = rev(c('BvsA', 'CvsA', 'BvsC')))
	g3 = ggplot(data, aes(x = ID, y = NES, group = comparison, fill = comparison)) + 
             geom_bar(stat = "identity", position=position_dodge(0.75)) + 
             coord_flip() +
	     scale_x_discrete(label = substring(levels(data$ID), 0, 60)) +
             theme_cowplot()
	pdf(paste0("./results/02-10-33-gsea-gobp-barplot-", 
	           names(gobp.gsea.dflist)[j], ".pdf"), 
	    width = 15, height = 10)
	print(g3)
	dev.off()
}

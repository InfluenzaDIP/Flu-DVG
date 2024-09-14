library(Seurat)
library(patchwork)
library(tidyverse)
library(magrittr)
library(future)

set.seed(1515)
options(future.globals.maxSize = 80000 * 1024^2, seed = TRUE)
plan("multicore", workers = 16)
plan()


# 1. Load count matrix from CellRanger
# step1 list sample directories ----------------------------------------------
sample <- list.dirs(path = '/workdir/licd/Jida/1-count/',
                    full.names = F,
                    recursive = F)
sample = sample[-4]
dir.ls = paste0( '/workdir/licd/Jida/1-count/', sample, '/outs/filtered_feature_bc_matrix')
names(dir.ls) <- sample
# step2 check whether dir exist -------------------------------------------
dir.ls %>% map( ~ dir.exists(.x))
# step3 create seurat per samples -----------------------------------------
obj.ls <- dir.ls %>% map( ~ Read10X(.x)) %>% map( ~ CreateSeuratObject(.x, min.cells = 0, min.features = 0))


# 2.
for(i in 1:length(obj.ls)){
  #给细胞barcode加个前缀，防止合并后barcode重名
  obj.ls[[i]] <- RenameCells(obj.ls[[i]], add.cell.id = sample[i])
  #加信息
  obj.ls[[i]][["sample"]]  = sample[i]
  #计算线粒体基因比例
  if(T){    
    obj.ls[[i]][["percent.mt"]] <- PercentageFeatureSet(obj.ls[[i]], pattern = "^mt-") 
  }
  #计算核糖体基因比例
  if(T){
    obj.ls[[i]][["percent.rb"]] <- PercentageFeatureSet(obj.ls[[i]], pattern = "^Rp[sl]")
  }
}
saveRDS(obj.ls, file = "/workdir/wangph/proj/20231102_Jida/SC-3sample/1.1.1.list-raw-20231102.rds")


# 3. Merge individuals
combined <- merge(x = obj.ls[[1]], y = obj.ls[2:length(obj.ls)])
combined
head(combined@meta.data)


# 4. virus
virus_B = read.csv("/workdir/wangph/proj/20231102_Jida/blast/B.output", sep = '\t')
virus_B = virus_B$Cell.barcode
virus_B = unique(virus_B)
virus_B = paste0('B_', virus_B, '-1')
head(virus_B)
virus_C = read.csv("/workdir/wangph/proj/20231102_Jida/blast/C.output", sep = '\t')
virus_C = virus_C$Cell.barcode
virus_C = unique(virus_C)
virus_C = paste0('C_', virus_C, '-1')
head(virus_C)
virus_D = read.csv("/workdir/wangph/proj/20231102_Jida/blast/D.output", sep = '\t')
virus_D = virus_D$Cell.barcode
virus_D = unique(virus_D)
virus_D = paste0('D_', virus_D, '-1')
head(virus_D)
virus = c(virus_B, virus_C, virus_D)

combined@meta.data$Virus = 'No'
combined@meta.data$Virus[rownames(combined@meta.data) %in% virus] <- 'Yes'
head(combined@meta.data)
table(combined@meta.data$Virus)
table(combined@meta.data$Virus, combined@meta.data$sample)


# 5. save
saveRDS(combined, file = "/workdir/wangph/proj/20231102_Jida/SC-3sample/1.1.2.merged-raw-20231102.rds")

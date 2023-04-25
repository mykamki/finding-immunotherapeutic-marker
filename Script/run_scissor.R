#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(Scissor)
library(data.table)
library(dplyr)
library(purrr)
indir_bulk <- args[1]
indir_sc <- args[2]
outdir <- args[3]


### 01. Load bulk data
load(file = paste0(blocation, "bulk_dataset.RData"))
load(file = paste0(blocation, "phenotype.RData"))

# prepare bulk dataset
bulk_name <- rownames(bulk_dataset)
bulk_dataset <- sapply(bulk_dataset, as.double)
rownames(bulk_dataset) <- bulk_name

# remove one patients
ix <- which(colnames(bulk_dataset) %in% "BACI165_2")
bulk_dataset <- bulk_dataset[,-ix]
phenotype <- phenotype[-ix,]


#### 2. Load scdata ####
#sclocation <- args[2]
load(file = args[2])
if ("data1" %in% ls()) {

  selection.method = "vst"
  resolution = 0.6
  dims_Neighbors = 1:10
  dims_TSNE = 1:10
  dims_UMAP = 1:10
  verbose = TRUE
  DefaultAssay(data1) <- "RNA"
  data1 <- FindVariableFeatures(object = data1, selection.method = selection.method, verbose = verbose)
  data1 <- ScaleData(object = data1, verbose = verbose)
  data1 <- RunPCA(object = data1, features = VariableFeatures(data1), verbose = verbose)
  data1 <- FindNeighbors(object = data1, dims = dims_Neighbors, verbose = verbose)
  data1 <- FindClusters( object = data1, resolution = resolution, verbose = verbose)
  data1 <- RunTSNE(object = data1, dims = dims_TSNE)
  data1 <- RunUMAP(object = data1, dims = dims_UMAP, verbose = verbose)
  save(data1, file = paste0(args[4], 'preprocessed_data1.RData'))
  sc_dataset <- data1
  
  } else {
  
  sc_dataset <- Seurat_preprocessing(count.data, verbose = F)
  save(sc_dataset, file = paste0(args[4], 'preprocessed_data1.RData'))
  }


sclocation <- paste0(paste(rev(rev(strsplit(args[2], "/")[[1]])[-1]), collapse = "/"), "/")
print(sclocation)

#### 3. Execute Scissor to select the informative cells ####
print(args[3])
alp <- ""
if (args[3] == "null") {
	alp <- NULL
} else {
	alp <- as.double(args[3])
}



print(alp)
infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = alp, 
                 family = "cox", Save_file = paste0(args[4], 'Scissor_Bladder_survival.RData'))

save(infos1, file = paste0(args[4], 'infos1.RData'))

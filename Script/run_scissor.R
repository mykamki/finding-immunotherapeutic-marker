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
load(file = paste0(indir_bulk, "bulk_dataset.RData"))
load(file = paste0(indir_bulk, "bulk_clinical.RData"))




### 02. Make phenotype data
phenotype <- bulk_clinical[,c("Patient", "overall survival", "alive")]
colnames(phenotype) <- c("ID","time","status")
phenotype$status <- ifelse(phenotype$status == "No", 1, 0)
phenotype$time <- as.double(phenotype$time)

# check whether bulk_dataset colname matched with phenotype ID
if (identical(phenotype$ID , colnames(bulk_dataset)))
  bulk_dataset <- bulk_dataset[,phenotype$ID]
}
rownames(phenotype) <- phenotype$ID
phenotype <- phenotype[,c("time","status")]



### 03. Load scdata 
load(file = paste0(indir_sc, "preprocessed_data1.RData"))



### 03. Preprocessing for scissor
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
sc_dataset <- data1


### 03. Execute Scissor to select the informative cells 
bulk_dataset <- as.matrix(bulk_dataset)
alp <- NULL
infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = alp, 
                 family = "cox", Save_file = paste0(outdir, 'Scissor_Bladder_survival.RData'))



### 04. Save the results
save(sc_dataset, file = paste0(outdir, 'sc_dataset.RData'))
save(infos1, file = paste0(outdir, 'infos1.RData'))



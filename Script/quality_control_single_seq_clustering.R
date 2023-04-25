#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(cluster, quietly = TRUE)
library(dplyr)
indir <- args[1]



#### 01. Load data ####
load(file = paste0(args[1], "data1.RData"))


#### 02. Cluster the cells ####
max.dim <- as.double(args[2])
opt.resol <- as.double(args[3])
data1 <- FindNeighbors(data1, dims = 1:max.dim) # default dims_Neighbors = 1:10
data1 <- FindClusters(data1, resolution = opt.resol)


#### 03. Run non-linear dimensional reduction (UMAP/tSNE) ####
# Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. 
# The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. 
# Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots. 
# As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.
data1 <- RunTSNE(object = data1, dims = 1:max.dim)
data1 <- RunUMAP(object = data1, dims = 1:max.dim)


save(data1, file = paste0(args[4], "preprocessed_data1.RData"))


#### 04. Plot ####
set.seed(12345)

pdf(file = paste0(args[4], "umap.pdf"))
DimPlot(data1, reduction = "umap")
dev.off()

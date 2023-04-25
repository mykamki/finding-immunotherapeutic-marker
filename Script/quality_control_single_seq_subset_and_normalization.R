#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
indir <- args[1]
outdir <- args[2]
min.nfeature <- as.double(args[3])
max.nfeature <- as.double(args[4])
max.mt <- as.double(args[5])



### 01. Load data 
load(file = paste0(indir, "sc_dataset.RData"))



### 02. Subset data
data1 <- subset(data, subset = nFeature_RNA > min.nfeature & nFeature_RNA < max.nfeature & percent.mt < max.mt)  # min.cells = 400, min.features = 0
cat("## The number of cell before QC : ", length(data$orig.ident), " ##", "\n")
cat("## The number of cell after QC : ", length(data1$orig.ident), " ##", "\n")



### 03. Normalizing the data 
data1 <- NormalizeData(data1, normalization.method = "LogNormalize", scale.factor = 10000)
data1 <- FindVariableFeatures(data1, selection.method = "vst", nfeatures = 2000)



### 04. Scaling the data 
all.genes <- rownames(data1)
data1 <- ScaleData(data1, features = all.genes)



### 05. Perform linear dimensional reduction
data1 <- RunPCA(data1, features = VariableFeatures(object = data1))



### 06. Save the results and elbow plot
save(data1, file = paste0(outdir, "normalized_data1.RData"))
pdf(file = paste0(outdir, "elbowplot.pdf"))
ElbowPlot(data1, ndims= 50)
dev.off()



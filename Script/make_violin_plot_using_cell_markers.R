#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(purrr)
library(gridExtra)
indir <- args[1]



### 01. load data
load(paste0(indir, "preprocessed_data1.RData"))
if ("sc_dataset" %in% ls()) {
  data1 <- sc_dataset
  DefaultAssay(object =data1) <- "RNA"
}

DefaultAssay(object =data1) <- "RNA"



### 02. make cell type result directory
ls <- system(paste0("ls ", indir), intern = T)
if ("cell_type_ann_result" %in% ls) {
  print("already result's directory existed")
} else {
  system(paste0("mkdir ",indir, "cell_type_ann_result"))
}
outdir <- paste0(indir, "cell_type_ann_result/")



### 03. Violin plot
# urothelial cells
png(file = paste0(outdir, "epithelialcell.png"), width=1200, height=1000)
p1 <- VlnPlot(data1, features = "KRT19", pt.size = 0)
p2 <- VlnPlot(data1, features = "KRT18", pt.size = 0)
gridExtra::grid.arrange(p1, p2, nrow = 2)
dev.off()


# myeloid/macrophage
png(file = paste0(outdir, "myeloid.macrophage.png"), width=1200, height=1000)
p1 <- VlnPlot(data1, features = "LYZ", pt.size = 0)
p2 <- VlnPlot(data1, features = "CSF1R", pt.size = 0)
gridExtra::grid.arrange(p1,p2, nrow = 2)
dev.off()


# T cell
png(file = paste0(outdir, "Tcell.png"), width=1200, height=1000)
p1 <- VlnPlot(data1, features = "CD3D", pt.size = 0)
p2 <- VlnPlot(data1, features = "CD3E", pt.size = 0)
gridExtra::grid.arrange(p1, p2, nrow = 2)
dev.off()


# endothelial cells
png(file = paste0(outdir, "endothelialcell.png"), width=1200, height=1000)
p1 <- VlnPlot(data1, features = "VWF", pt.size = 0)
gridExtra::grid.arrange(p1, nrow = 2)
dev.off()


# B cell
png(file = paste0(outdir, "Bcell.png"), width=1200, height=1000)
p1 <- VlnPlot(data1, features = "CD79A", pt.size = 0)
p2 <- VlnPlot(data1, features = "MZB1", pt.size = 0)
gridExtra::grid.arrange(p1, p2, nrow = 2)
dev.off()


# fibroblasts
png(file = paste0(outdir, "fibroblast.png"), width=1200, height=1000)
p1 <- VlnPlot(data1, features = "TAGLN", pt.size = 0)
p2 <- VlnPlot(data1, features = "PDPN", pt.size = 0)
gridExtra::grid.arrange(p1, p2, nrow = 2)
dev.off()

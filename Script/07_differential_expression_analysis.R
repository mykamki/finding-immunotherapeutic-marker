#!/usr/bin/env Rscript
# 07_differential_expression_analysis.R $INDIR 


### 01. Setting
source("Script/Functions/func_exploring.R")
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(dplyr)
library(DESeq2)


### 02. Load data
load(paste0(indir, "infos1.RData"))
load(paste0(indir, "annotated_data1.RData"))



### 03. Extract Scissor cells
pos_cell_names <- infos1$Scissor_pos
neg_cell_names <- infos1$Scissor_neg



### 04. Identify scissor cells 
celltypes <- as.vector(Idents(data1))
celltypes[names(Idents(data1)) %in% neg_cell_names] <- "neg_scissor"
celltypes[names(Idents(data1)) %in% pos_cell_names ] <- "pos_scissor"
Idents(data1) <- celltypes



### 05. Do DE analysis
myident1 <- "neg_scissor"
others <- c("T_cell", "B_cell", "Myeloid_macrophage", "Fibroblast", "other", "Endothelial_cell")  
de.neg <- FindMarkers(data1, ident.1 = myident1, ident.2 = others, test.use = "DESeq2")
de.neg <- test_de_cutoff(de.neg, 1.5)
save(de.neg, file = paste0(indir, "negcell.de.result.RData"))
write.csv(de.neg, file = paste0(indir, "scissor_negative_de_result.csv"), row.names = F, quote = F)


### 06. Vocano plot

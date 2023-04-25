#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(data.table)
library(dplyr)
indir <- args[1]
cell_ann <- args[2]
targetcell <- args[3]
outdir <- args[4]



### 01. Load data
load(paste0(indir, "preprocessed_data1.RData"))
source(cell_ann)



### 02. Change cell identity classes
cell.name <- names(my_cluster_ann)
data1.new.idents <- case_when(
                      Idents(data1) %in% my_cluster_ann[[1]] ~  cell.name[1],
                      Idents(data1) %in% my_cluster_ann[[2]] ~  cell.name[2],
                      Idents(data1) %in% my_cluster_ann[[3]] ~  cell.name[3],
                      Idents(data1) %in% my_cluster_ann[[4]] ~  cell.name[4],
                      Idents(data1) %in% my_cluster_ann[[5]] ~  cell.name[5],
                      Idents(data1) %in% my_cluster_ann[[6]] ~  cell.name[6],
                      TRUE ~ NA_character_
                      )
Idents(data1) <- data1.new.idents



### 03. Subset cells 
cat("####################################################\n")
cat("#### Before removing ", targetcell, " : ", length(data1$orig.ident), "####\n")
data1 <- subset(data1, idents = targetcell, invert = TRUE)

cat("##### After removing ", targetcell, " : ", length(data1$orig.ident), "#####\n")
cat("####################################################\n")



### 04. Save the results
save(data1, file = paste0(outdir, "preprocessed_data1.RData"))



library(Seurat)
library(dplyr)
library(DESeq2)


### 
 Idents(data1)

myident1 <- "neg_scissor"
others <- c("T_cell", "B_cell", "Myeloid_macrophage", "Fibroblast", "other", "Endothelial_cell")  
de.neg <- FindMarkers(data1, ident.1 = myident1, ident.2 = others, test.use = "DESeq2")
de.neg <- test_de_cutoff(de.neg, 1.5)

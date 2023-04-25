#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(data.table)
library(dplyr)
indir <- args[1]
outdir <- args[2]


### 01. Load bulk data
bulk_dataset <- fread(paste0(indir, "GSE176307_baci_rsem_RS_BACI_headers_tab.txt.gz"), header = T, data.table = F)
bulk_map <- fread(paste0(indir, "GSE176307_BACI_Omniseq_Sample_Name_Key_submitted_GEO_v2.csv.gz"), data.table = F)
bulk_dataset <- bulk_dataset[-1,]



### 02. Match column ID 
colnames(bulk_map) <- c("id", "rnaid")
bulk_map[47,1] <- "BACI165_1"
bulk_map[47,2] <- "RS-03239001"
bulk_map[90,1] <- "BACI165_2"
bulk_map[90,2] <- "RS-03238964"

# sort patients
a <- data.frame(colnames(bulk_dataset)[2:ncol(bulk_dataset)], seq(1,90))
colnames(a) <- c("rnaid", "no")
a <- merge(a, bulk_map, by = "rnaid")
a <- a %>% arrange(no)
colnames(bulk_dataset) <- c("V1", a$id)
rownames(bulk_dataset) <- bulk_dataset$V1
bulk_dataset <- bulk_dataset[,-1]



### 03. Remove duplicated patients
ix <- which(colnames(bulk_dataset) %in% "BACI165_2")
bulk_dataset <- bulk_dataset[,-ix]



### 04. Save the results
save(bulk_dataset, file = paste0(outdir, "bulk_dataset.RData"))



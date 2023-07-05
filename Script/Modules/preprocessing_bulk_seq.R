#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(data.table)
library(dplyr)
source("Script/Functions/func_preprocessing.R")
indir <- args[1]
outdir <- args[2]


### 01. Load bulk data
bulk_metadata <- fread(paste0(indir, "GSE176307_series_matrix.txt.gz"), sep = "\t", fill = TRUE) 
bulk_dataset <- fread(paste0(indir, "GSE176307_BACI_tpm_gene.matrix.tsv.gz"))



bulk_dataset <- fread(paste0(indir, "GSE176307_baci_rsem_RS_BACI_headers_tab.txt.gz"), header = T, data.table = F)
bulk_map <- fread(paste0(indir, "GSE176307_BACI_Omniseq_Sample_Name_Key_submitted_GEO_v2.csv.gz"), data.table = F)
bulk_dataset <- bulk_dataset[-1,]





### 03. Preprocessing clinical data
# The list of values for extract the metadataset
values <- c("Patient", "gender", "race", "ethnicity", "age", "smoking status" , "histology", "histology.other.description", 
  "t.stage.at.diagnosis", "n.stage.at.diagnosis" , "d.stage.at.diagnosis", "pack.years" , "overall survival",
  "alive", "primary tumor location", "io.lot" , "ecog", "met.lung.cumul", "met.liver.cumul", "met.bone.cumul",
  "met.other.cumul", "fgfr.mutation.yn" , "io.therapy", "duration.of.io.tx" , "io.response", "pfs" , "progressed",
   "tmb", "tmb.interpretation")

# extract clinical data
bulk_clinical <- as.data.frame(bulk_metadata[bulk_metadata$V1 %in% c("!Sample_title" , "!Sample_characteristics_ch1" ),])
bulk_clinical <- sapply(values, extract_value2list)
bulk_clinical <- t(do.call(rbind.data.frame, bulk_clinical))
colnames(bulk_clinical) <- values
rownames(bulk_clinical) <- as.vector(bulk_clinical[,1])
bulk_clinical <- as.data.frame(bulk_clinical)

# sort patients
a <- data.frame(colnames(bulk_dataset), seq(1,90))
colnames(a) <- c("Patient", "no")
a <- merge(a, bulk_clinical, by = "Patient")
bulk_clinical <- a %>% arrange(no)




### 04. Make phenotype data for scissor 
bulk_survival <- bulk_clinical[, c("Patient", "overall survival", "alive")]
colnames(bulk_survival) <- c("id", "time", "status")
bulk_survival$time <- as.double(bulk_survival$time)
bulk_survival$status <- ifelse(bulk_survival$status == "No", 1, 0)



### 05. Remove duplicated patients
ix <- which(colnames(bulk_dataset) %in% "BACI165_2")
bulk_dataset <- bulk_dataset[,-ix]
phenotype <- phenotype[-ix,]



### 06. Save data ####
save(bulk_dataset, file = paste0(outdir, "bulk_dataset.RData"))
save(phenotype, file = paste0(args[2], "phenotype.RData"))
save(bulk_clinical, file = paste0(args[2], "bulk_clinical.RData"))



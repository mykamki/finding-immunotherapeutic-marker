#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(data.table)
library(dplyr)
source("Script/Functions/func_preprocessing.R")
indir <- args[1]
outdir <- args[2]


### 01. Load bulk data
bulk_metadata <- fread(paste0(indir, "GSE176307_series_matrix.txt.gz"), sep = "\t", fill = TRUE) 
bulk_dataset <- fread(paste0(indir, "GSE176307_BACI_tpm_gene.matrix.tsv.gz"), data.table = F)



### 02. Preprocessing clinical data
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

# remove dulplicate patient (BACI165)
#grep("165",bulk_clinical$Patient) #47 48
bulk_clinical <- bulk_clinical[-48,]
bulk_clinical$Patient <- ifelse(bulk_clinical$Patient == "BACI165_1", "BACI165", bulk_clinical$Patient)



### 03. Preprocessing tpm data
# sort patients by clinical data order
genenames <- bulk_dataset$V1
bulk_dataset <- bulk_dataset[,bulk_clinical$Patient]
rownames(bulk_dataset) <- genenames



### 04. Save data ####
save(bulk_dataset, file = paste0(outdir, "bulk_dataset.RData"))
save(bulk_clinical, file = paste0(outdir, "bulk_clinical.RData"))



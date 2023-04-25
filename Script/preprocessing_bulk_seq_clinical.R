#### Prepare the bulk phenotype ####
library(data.table)
library(dplyr)
source("Functions/func_preprocessing.R")
indir <- "../Data/"
outdir <- "../Output/"



### 01. Load bulk data
bulk_metadata <- fread(paste0(indir, "GSE176307_series_matrix.txt.gz"), sep = "\t", fill = TRUE) 



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

# sort patients
a <- data.frame(colnames(bulk_dataset), seq(1,90))
colnames(a) <- c("Patient", "no")
a <- merge(a, bulk_clinical, by = "Patient")
bulk_clinical <- a %>% arrange(no)



### 03. Make phenotype data for scissor 
bulk_survival <- bulk_clinical[, c("Patient", "overall survival", "alive")]
colnames(bulk_survival) <- c("id", "time", "status")
bulk_survival$time <- as.double(bulk_survival$time)
bulk_survival$status <- ifelse(bulk_survival$status == "No", 1, 0)



### 04. Save data ####
save(phenotype, file = paste0(args[2], "phenotype.RData"))
save(bulk_clinical, file = paste0(args[2], "bulk_clinical.RData"))



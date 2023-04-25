#### Prepare the bulk phenotype ####
library(data.table)
library(dplyr)
indir <- ""
outdir <- ""



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


a <- data.frame(colnames(bulk_dataset), seq(1,90))
colnames(a) <- c("Patient", "no")
a <- merge(a, bulk_clinical, by = "Patient")
bulk_clinical <- a %>% arrange(no)

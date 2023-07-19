### 01. Setiing
library(limma)
library(dplyr)
indir <- args[1]
outdir <- args[2]


### 02. Load data
load(file = paste0(indir, "bulk_dataset.RData"))
load(file = paste0(indir, "bulk_clinical.RData"))


### 02. Make clinical_gse176307 data
clinical_gse176307 <- bulk_clinical[,c("Patient", "overall survival", "alive", "tmb", "io.response","gender", "race", "smoking status", "ecog")]
colnames(clinical_gse176307) <- c("ID", "time", "status", "tmb", "response", "sex", "ethnicity", "smk", "ecog")

clinical_gse176307$status <- ifelse(clinical_gse176307$status == "No", 1, 0)
clinical_gse176307$time <- as.double(clinical_gse176307$time)
clinical_gse176307$tmb <- as.double(clinical_gse176307$tmb)
clinical_gse176307$response <- factor(clinical_gse176307$response, levels = c("CR", "PR", "SD", "PD", "NA"))
clinical_gse176307$sex <- ifelse(clinical_gse176307$sex == "M", 1, 2)
clinical_gse176307$ethnicity <- ifelse(clinical_gse176307$ethnicity == "White", 1, 
						                      ifelse(clinical_gse176307$ethnicity == "Black or African American", 2, 
								                    ifelse(clinical_gse176307$ethnicity == "Asian", 3,
										                  ifelse(clinical_gse176307$ethnicity == "Unknown", 5, 4))))
clinical_gse176307$smk <- ifelse(clinical_gse176307$smk == "Current", 1, 
						                ifelse(clinical_gse176307$smk == "Light", 2, 
								              ifelse(clinical_gse176307$smk == "Former", 3,
										            ifelse(clinical_gse176307$smk == "Never", 4, NA))))
clinical_gse176307$ecog <- ifelse(clinical_gse176307$ecog == "0", 0, 
						                ifelse(clinical_gse176307$ecog == "1", 1, 
								              ifelse(clinical_gse176307$ecog == "2", 2,
										            ifelse(clinical_gse176307$ecog == "3", 3, NA))))

#sex # Male 1 Female 2
#ethnicity # white 1 Black or African American 2 Asian 3 Others 4 Unknown 5
#smk # Current smoker 1 Light smoker 2 Former smoker 3 Never smoker 4
#ecog # 0 1 2 3 NA


### 03. preprocessing bulk dataset for bladder signature
### ...03-1. make normalized zscore
#mygene <- c("SSR4", "CD74", "HLA-DPA1", "JCHAIN", "HLA-DRA", "RGS1")
mygene <- c("SSR4", "RGS1" ,"HLA-DRB5", "APOE", "C1QB",  "C1QA",  "APOC1",  "JCHAIN",  "C1QC", "DERL3")
v <- voom(bulk_dataset)
log_dataset <- v$E
mygene <- ifelse(mygene %in% "JCHAIN", "IGJ", mygene)
log_dataset_gse176307 <- apply(log_dataset[rownames(log_dataset) %in% mygene ,], 1, zscore_transform)
log_dataset_gse176307 <- t(log_dataset_gse176307)


### ...03-2. divide group by bladder signature
if (rownames(log_dataset_gse176307) %in% mygene %>% sum() == 10) {
	clinical_gse176307 <- make_genesignature(clinical_gse176307, log_dataset_gse176307)
}



### 04. Save data
save(clinical_gse176307, file = paste0(outdir, "clinical_gse176307.RData"))
save(log_dataset_gse176307, file = paste0(outdir, "log_dataset_gse176307.RData"))


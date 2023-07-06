### 01. Setiing
library(limma)
library(dplyr)


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
mygene <- c("SSR4", "CD74", "HLA-DPA1", "JCHAIN", "HLA-DRA", "RGS1")
zscore_transform <- function(x) {
  (x-mean(x))/sd(x)
}
v <- voom(bulk_dataset)
log_dataset <- v$E
log_dataset2 <- apply(log_dataset, 2, zscore_transform)



# check whether bulk_dataset colname matched with phenotype ID
if (identical(phenotype$ID , colnames(bulk_dataset))) {
  bulk_dataset <- bulk_dataset[,phenotype$ID]
}
rownames(phenotype) <- phenotype$ID
phenotype <- phenotype[,c("time","status")]


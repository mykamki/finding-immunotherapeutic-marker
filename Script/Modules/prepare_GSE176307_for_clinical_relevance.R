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
zscore_transform <- function(x) {
  (x-mean(x))/sd(x)
}
v <- voom(bulk_dataset)
log_dataset <- v$E
log_dataset2 <- apply(log_dataset, 2, zscore_transform)


### ...03-2. divide group by bladder signature
mygene <- ifelse(mygene %in% "JCHAIN", "IGJ", mygene)
rownames(log_dataset2) %in% mygene %>% sum() # check gene names
res2 <- apply(log_dataset2[rownames(log_dataset2) %in% mygene ,],2,mean)
sumgsig <- summary(res2)
names(res2)[which(res2<=sumgsig[3])] -> low # low ID
names(res2)[which(res2>sumgsig[3])] -> high # high ID

clinical_gse176307 <- clinical_gse176307 %>% mutate(Novel_Signature = ifelse(ID %in% high, "High", "Low"))
clinical_gse176307$Novel_Signature <- factor(clinical_gse176307$Novel_Signature)

if (identical(clinical_gse176307$ID, names(res2))) {
	clinical_gse176307$Novel_Signature_score <- res2
} else {
	res2 <- res2[clinical_gse176307$ID]
	clinical_gse176307$Novel_Signature_score <- res2
}



### 04. Save data
log_dataset_gse176307 <- log_dataset2
save(clinical_gse176307, file = paste0(outdir, "clinical_gse176307.RData"))
save(log_dataset_gse176307, file = paste0(outdir, "log_dataset_gse176307.RData"))



### 01. Setting
# conda activate imvcore
library(IMvigor210CoreBiologies)
library(limma)
library(dplyr)
source("/Script/Functions/func_clinical_relevance.R")



### 02. Preprocessing bulk_dataset
### ... 02-1. change count to tpm
data(cds) # load data
imvCounts <- as.data.frame(counts(cds)) %>% as.matrix() # count
imvTpms <- apply(imvCounts, 2, function(x) tpm(x, fData(cds)$length)) # tpm

### ... 02-2. change tpm to normalized zscore
v <- voom(imvTpms)
log_dataset <- v$E
log_dataset_imvigor210core <- apply(log_dataset, 2, zscore_transform)


### 03. Make clinical data
#sex # Male 1 Female 2
#ethnicity # white 1 Black or African American 2 Asian 3 Others 4 Unknown 5
#smk # Current smoker 1 Light smoker 2 Former smoker 3 Never smoker 4
#ecog # 0 1 2 3 NA

pData(cds) -> pdata_imvigor210core
clinical_imvigor210core <- pdata_imvigor210core[,c("FMOne mutation burden per MB","Best Confirmed Overall Response", 
                                                   "Sex", "Race", "Tobacco Use History", "Baseline ECOG Score", 
                                                   "censOS", "os")]
colnames(clinical_imvigor210core) <- c("tmb", "response", "sex", "ethnicity", "smk", "ecog","status", "time")
clinical_imvigor210core$ID <- rownames(pdata_imvigor210core)
clinical_imvigor210core$sex <- ifelse(clinical_imvigor210core$sex == "M", 1, 2)
clinical_imvigor210core$ethnicity <- ifelse(clinical_imvigor210core$ethnicity == "WHITE", 1, 
						ifelse(clinical_imvigor210core$ethnicity == "BLACK OR AFRICAN AMERICAN", 2, 
								ifelse(clinical_imvigor210core$ethnicity == "ASIAN", 3,
										ifelse(clinical_imvigor210core$ethnicity == "UNKNOWN", 5, 4))))
clinical_imvigor210core$smk <- ifelse(clinical_imvigor210core$smk == "CURRENT", 1, 
						ifelse(clinical_imvigor210core$smk == "Light", 2, 
								ifelse(clinical_imvigor210core$smk == "PREVIOUS", 3,
										ifelse(clinical_imvigor210core$smk == "NEVER", 4, NA))))
clinical_imvigor210core$ecog <- ifelse(clinical_imvigor210core$ecog == "0", 0, 
						ifelse(clinical_imvigor210core$ecog == "1", 1, 
								ifelse(clinical_imvigor210core$ecog == "2", 2,
										ifelse(clinical_imvigor210core$ecog == "3", 3, NA))))
clinical_imvigor210core$time <- clinical_imvigor210core$time *30.436875          

		 

### 04. preprocessing bulk dataset for bladder signature
### ...04-1. make normalized zscore
#mygene <- c("SSR4", "CD74", "HLA-DPA1", "JCHAIN", "HLA-DRA", "RGS1")
mygene <- c("SSR4", "RGS1" ,"HLA-DRB5", "APOE", "C1QB",  "C1QA",  "APOC1",  "JCHAIN",  "C1QC", "DERL3")                 
#fData(cds)[fData(cds)$Symbol %in% mygene,]$entrez_id
                 
res2 <- apply(log_dataset_imvigor210core[rownames(log_dataset_imvigor210core) %in% fData(cds)[fData(cds)$Symbol %in% mygene,]$entrez_id,],2,mean)
sumgsig <- summary(res2)
names(res2)[which(res2<=sumgsig[3])] -> low # low ID
names(res2)[which(res2>sumgsig[3])] -> high # high ID

### ... 04-2. Divide group 
clinical_imvigor210core <- clinical_imvigor210core %>% mutate(Novel_Signature = ifelse(ID %in% high, "High", "Low"))
clinical_imvigor210core$Novel_Signature <- factor(clinical_imvigor210core$Novel_Signature)

if (identical(clinical_imvigor210core$ID, names(res2))) {
	res2 <- res2[clinical_imvigor210core$ID]
	}
clinical_imvigor210core$Novel_Signature_score <- res2

		 

### 04. Save data
save(clinical_imvigor210core, file = paste0(outdir, "clinical_imvigor210core.RData"))
save(log_dataset_imvigor210core, file = paste0(outdir, "log_dataset_imvigor210core.RData"))

		 

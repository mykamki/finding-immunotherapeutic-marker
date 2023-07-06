
### 01. Setting
# conda activate imvcore
library(IMvigor210CoreBiologies)
library(limma)
library(dplyr)


### 02. Preprocessing bulk_dataset
# make functions

# preprocessing rna-seq dataset
# count
data(cds)
imvCounts <- as.data.frame(counts(cds)) %>% as.matrix()

# tpm
imvTpms <- apply(imvCounts, 2, function(x) tpm(x, fData(cds)$length))

# zscore
v <- voom(imvTpms)
log_dataset <- v$E
log_dataset_imvigor210core <- apply(log_dataset, 2, zscore_transform)

### make gene group
mygene <- c("SSR4", "CD74", "HLA-DPA1", "JCHAIN", "HLA-DRA", "RGS1", "IGJ")
#fData(cds)[fData(cds)$Symbol %in% mygene,]$entrez_id

pData(cds) -> clinical_imvigor210core

clinical_imvigor210core$ID <- rownames(clinical_imvigor210core)
clinical_imvigor210core$time <- clinical_imvigor210core$os * 30.436875
clinical_imvigor210core$status <- clinical_imvigor210core$censOS
clinical_imvigor210core$response <- clinical_imvigor210core$`Best Confirmed Overall Response`
clinical_imvigor210core$tmb <- clinical_imvigor210core$`FMOne mutation burden per MB`

res2 <- apply(log_dataset_imvigor210core[rownames(log_dataset_imvigor210core) %in% fData(cds)[fData(cds)$Symbol %in% mygene,]$entrez_id,],2,mean)
sumgsig <- summary(res2)
names(res2)[which(res2<=sumgsig[3])] -> low # low ID
names(res2)[which(res2>sumgsig[3])] -> high # high ID
clinical_imvigor210core <- clinical_imvigor210core %>% mutate(Gene_group = ifelse(ID %in% high, "High", "Low"))
identical(clinical_imvigor210core$ID, names(res2)) # TRUE
clinical_imvigor210core$genesignature <- res2

### save
save(clinical_imvigor210core, file = "clinical_imvigor210core.RData")
save(log_dataset_imvigor210core, file = "log_dataset_imvigor210core.RData")

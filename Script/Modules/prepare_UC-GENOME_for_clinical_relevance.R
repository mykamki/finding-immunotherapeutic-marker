### 01. Setting
library(cgdsr)
library(dplyr)
source("/Script/Functions/func_clinical_relevance.R")
outdir <- args[1]



### 02. Make UC-GENOME object from cgdsr
###...02-1. Loading data (UC-GENOME package) ###
mycgds <- CGDS("http://www.cbioportal.org/") # make cgds object

###...02-2. Get target study id 
cancerstudy <- getCancerStudies(mycgds) # get cencer studies
idx <- grep("Urothelial", cancerstudy$name, ignore.case= T) # find ovarian cancer study
cancerstudy$name[idx] 
mycancerstudy <- cancerstudy[idx[11],1] #get Urothelial cancer study id

###...02-3. Extract samples id
bladdercase <- getCaseLists(mycgds, mycancerstudy) # get case info object of Urothelial cancer data
bladdercase[,2]
mycaselist <- bladdercase[2,1] # get "Samples with mRNA data (RNA Seq V2)"

###...02-4. Extract features id
myfeature <- getGeneticProfiles(mycgds, mycancerstudy)
myfeature[,2]
myprofile <- myfeature[2,1] # "UC_GENOME_RNAseq_zscore"
                            # mRNA expression UQ_Log2_zscore

###...02-5. Save objects
uc_genome_object <- list(myprofile, mycaselist)
save(uc_genome_object, file = paste0(outdir, "uc_genome_object.RData"))



### 03. Make clinical data
#sex # Male 1 Female 2
#ethnicity # white 1 Black or African American 2 Asian 3 Others 4 Unknown 5
#smk # Current smoker 1 Light smoker 2 Former smoker 3 Never smoker 4
#ecog # 0 1 2 3 NA
pdata_ucgenome <- getClinicalData(mycgds, mycaselist)
                 
clinical_ucgenome <- pdata_ucgenome[,c("BEST_RESPONSE_IMMUNOTHERAPY", "SEX", "RACE", "SMOKING_STATUS", 
                                       "BASELINE_ECOG", "SURVIVAL_STATUS", "IMMUNO_SURVIVAL", "MUTATION_COUNT")]
colnames(clinical_ucgenome) <- c("response","sex", "ethnicity", "smk", "ecog","status", "time", "tmb")
clinical_ucgenome$ID <- rownames(pdata_ucgenome)
clinical_ucgenome$sex <- ifelse(clinical_ucgenome$sex == "Male", 1, 2)
clinical_ucgenome$ethnicity <- ifelse(clinical_ucgenome$ethnicity == "White", 1, 
						ifelse(clinical_ucgenome$ethnicity == "Black or African America", 2, 
								ifelse(clinical_ucgenome$ethnicity == "Asian", 3,
										ifelse(clinical_ucgenome$ethnicity == "Unknown", 5, 4))))
clinical_ucgenome$smk <- ifelse(clinical_ucgenome$smk == "Current Smoker", 1, 
						ifelse(clinical_ucgenome$smk == "Light", 2, 
								ifelse(clinical_ucgenome$smk == "Former Smoker", 3,
										ifelse(clinical_ucgenome$smk == "Never Smoker", 4, NA))))
clinical_ucgenome$ecog <- ifelse(clinical_ucgenome$ecog == "0", 0, 
						ifelse(clinical_ucgenome$ecog == "1", 1, 
								ifelse(clinical_ucgenome$ecog == "2", 2,
										ifelse(clinical_ucgenome$ecog == "3", 3, NA))))
clinical_ucgenome$status <- ifelse(clinical_ucgenome$status == "Dead", 1, 0)
clinical_ucgenome$time <- clinical_ucgenome$time *30.436875
clinical_ucgenome$response <- 
ifelse(clinical_ucgenome$response == "Partial Response", "PR",
				ifelse(clinical_ucgenome$response == "Not Evaluable", "NE",
				ifelse(clinical_ucgenome$response == "Complete Response", "CR",
				ifelse(clinical_ucgenome$response == "Progressive Disease", "PD",	
				ifelse(clinical_ucgenome$response == "Unknown", "NA",	
				ifelse(clinical_ucgenome$response == "Stable Disease", "SD","NA")
					)
				)
			)
		)
)
clinical_ucgenome$response <- factor(clinical_ucgenome$response ,
levels = c("CR", "PR", "SD", "PD" ,"NA","NE"))




### 04. preprocessing bulk dataset for bladder signature
### ...04-1. load mrna expression matrix
#mygene <- c("SSR4", "CD74", "HLA-DPA1", "JCHAIN", "HLA-DRA", "RGS1")
mygene <- c("SSR4", "RGS1" ,"HLA-DRB5", "APOE", "C1QB",  "C1QA",  "APOC1",  "JCHAIN",  "C1QC", "DERL3") 
log_dataset_ucgenome <- getProfileData(mycgds, mygene, myprofile, mycaselist)
log_dataset_ucgenome <- t(log_dataset_ucgenome )

### ...04-2. Check analysis target patients
apply(apply(log_dataset_ucgenome, 2, is.na),2, any) %>% sum() #the number who has NA in 10 genes at least one more
colnames(log_dataset_ucgenome)[apply(apply(log_dataset_ucgenome, 2, is.na),2, any)] -> na_pid # patient ID

pdata_ucgenome %>% 
filter(IMMUNOTHERAPY == "Yes") %>% nrow() # 125

pdata_ucgenome %>% 
filter(!rownames(pdata_ucgenome) %in% na_pid) %>% 
filter(IMMUNOTHERAPY == "Yes") %>% # 121
filter(SURVIVAL_STATUS %in% c("Alive", "Dead")) %>% # 108
filter(!is.na(IMMUNO_SURVIVAL)) %>% rownames() -> target_pid  # 101

# Total number of patients : 180
# The number of patients who were taken immnotherapy : 125
# The number of patients who were taken immnotherapy with all 10 gene expression data : 121
# The number of patients who were taken immnotherapy with all 10 gene expression data with survival data : 101

### ...04-3. clinical analysis data
clinical_ucgenome <- clinical_ucgenome[clinical_ucgenome$ID %in% target_pid,]

if (identical(clinical_ucgenome$ID, colnames(log_dataset_ucgenome))) {
  print("good")
  } else {
  log_dataset_ucgenome <- log_dataset_ucgenome[,clinical_ucgenome$ID]
  }
mygene <- gsub("-",".", mygene, fixed = T)

###... 04-4. Divide patients
if (rownames(log_dataset_ucgenome) %in% mygene %>% sum() == 10) {
	clinical_ucgenome <- make_genesignature(clinical_ucgenome, log_dataset_ucgenome)
}

log_dataset_ucgenome



		 
### 05. Save data
save(clinical_ucgenome, file = paste0(outdir, "clinical_ucgenome.RData"))
save(log_dataset_ucgenome, file = paste0(outdir, "log_dataset_ucgenome.RData"))

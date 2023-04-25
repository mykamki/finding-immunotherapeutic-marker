#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
indir <- args[1]
outdir <- args[2]



### 01. Load data
dic_list <- c("SC1", "SC2", "SC3", "SC4", "SC5", "SC6", "SC7")
for (n in 1:length(dic_list)) {
  print(paste0("sc",n)) 
  load(paste0(location, dic_list[n], "/preprocessed_data1.RData"))
  data1 <- SCTransform(data1, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
  assign(paste0("sc",n), data1)
}

remove(data1)



### 02. make list of sc dataset
sc_list <- list(sc1 = sc1,
                sc2 = sc2,
                sc3 = sc3,
                sc4 = sc4,
                sc5 = sc5,
                sc6 = sc6,
                sc7 = sc7)
  


### 03. Pre-processing data for integration
# Select the most variable features to use for integration
n.feature <- 3000              
integ_features <- SelectIntegrationFeatures(object.list = sc_list, 
                                            nfeatures = n.feature) 

# Prepare the SCT list object for integration
sc_list <- PrepSCTIntegration(object.list = sc_list, 
                                   anchor.features = integ_features)
                                   
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = sc_list, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
                                        
# Integrate across conditions
data1 <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")   



### 04. Save the results
save(data1, file = paste0(ourdir, "data1.RData"))



                                   

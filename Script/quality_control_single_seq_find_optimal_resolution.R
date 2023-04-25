#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
source("Functions/func_preprocessing.R")
library(cluster, quietly = TRUE)
library(dplyr)
library(boot)
library(ggplot2)
indir <- args[1]
max.dim <- as.double(args[2])
outdir <- indir



### 01. load data
load(file = paste0(indir, "normalized_data1.RData"))



### 02. Cluster the cells
data1 <- FindNeighbors(data1, dims = 1:max.dim) # default dims_Neighbors = 1:10



### 03. find resolution
# We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cell.
# Optimal resolution often increases for larger datasets.
resolutionrange <- c(0.5, 0.8, 1.0, 1.2, 1.6, 2.0, 4.0, 6.0, 8.0, 12.0, 16.0)
total_sil_df <- data.frame()
boot_res <- data.frame()

for (n in 1:length(resolutionrange)) {
  resol <- resolutionrange[n]
  data1 <- FindClusters(data1, resolution = resol)
  
  # calculate silhouette metric
  dist.matrix <- dist(x = Embeddings(object = data1[["pca"]])[, 1:max.dim])
  clusters <- data1$seurat_clusters
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
  
  sil_df <- data.frame(data1$seurat_clusters, sil[,3], rep(resol, length(sil[,3])))
  colnames(sil_df) <- c("clu","sil", "resol")
  cat("\n################################\n")
  cat("Finished! ", resol, " resolution calculating\n")
  
  # boot for confidence interval
  cluster.sil.summary <- sil_df %>% group_by(resol, clu) %>% summarize(n = mean(sil))
  boot_out <- boot(data = cluster.sil.summary$n, statistic = med, R = 2500)
  boot_ci <- boot.ci(boot_out, conf = 0.95, type = "bca")
  
  boot_res[n,1] <- resol
  boot_res[n,2] <- boot_ci$bca[4]
  boot_res[n,3] <- boot_ci$t0
  boot_res[n,4] <- boot_ci$bca[5]
  
  total_sil_df <- rbind(total_sil_df, sil_df)
  }

cat("################################\n")
colnames(boot_res) <- c("Resolution", "low_med", "med", "high_med")  
threshold <- max(boot_res$low_med)



### 04. Plot boxplot of Silhouette score per each resolution
cluster.sil.summary <- total_sil_df %>% group_by(resol, clu) %>% summarize(n = mean(sil))
p <- ggplot(cluster.sil.summary, aes(x=factor(resol), y=n, fill = factor(resol))) + 
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  scale_fill_brewer(palette="Set3") +
  theme_light() +
  xlab("Resolution") +
  ylab("Silhouette score") +
  geom_hline(yintercept= threshold, linetype="dashed", color = "red")



### 05. Save the results
save(boot_res, file = paste0(outdir, "boot_res.RData"))  
save(total_sil_df, file = paste0(outdir, "total_sil_df.RData"))

pdf(file = paste0(outdir, "silhouette.pdf"))
p
dev.off()




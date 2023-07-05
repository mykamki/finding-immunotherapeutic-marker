#!/usr/bin/env Rscript
# 07_differential_expression_analysis.R $INDIR 


### 01. Setting
source("Script/Functions/func_exploring.R")
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(grid)
library(gridExtra)



### 02. Load data
load(paste0(indir, "infos1.RData"))
load(paste0(indir, "annotated_data1.RData"))



### 03. Extract Scissor cells
pos_cell_names <- infos1$Scissor_pos
neg_cell_names <- infos1$Scissor_neg



### 04. Identify scissor cells 
celltypes <- as.vector(Idents(data1))
celltypes[names(Idents(data1)) %in% neg_cell_names] <- "neg_scissor"
celltypes[names(Idents(data1)) %in% pos_cell_names ] <- "pos_scissor"
Idents(data1) <- celltypes



### 05. Do DE analysis
myident1 <- "neg_scissor"
others <- c("T_cell", "B_cell", "Myeloid_macrophage", "Fibroblast", "other", "Endothelial_cell")  
de.neg <- FindMarkers(data1, ident.1 = myident1, ident.2 = others, test.use = "DESeq2")
de.neg <- test_de_cutoff(de.neg, 1.5)
save(de.neg, file = paste0(indir, "negcell.de.result.RData"))
write.csv(de.neg, file = paste0(indir, "scissor_negative_de_result.csv"), row.names = F, quote = F)



### 06. Volcano plot
volcano_p <- my_volcano(de.neg)

png(file = paste0(indir, "scissor_negative_de_vocanoplot.png"), width = 500, height= 500)
volcano_p
dev.off()


### 07. DE gene Violin plot
###...07-1. make de gene list 
de.neg %>% filter(expression == "Up-regulated") %>% rownames() -> degenes

# "SSR4" "CD74" "HLA-DPA1" "JCHAIN" "HLA-DRA"  "RGS1"
# "SSR4" "RGS1" "HLA-DRB5" "APOE" "C1QB"  "C1QA"  "APOC1"  "JCHAIN"  "C1QC"  "DERL3"   

degenes <- degenes[grep("IG", degenes, invert = T)]

###...07-2. make group (scissor neg cell vs. others) 
cell_group_df <- as.data.frame(Idents(data1))
colnames(cell_group_df) <- "celltype"
cell_group_df$group1 <- ifelse(cell_group_df$celltype == "neg_scissor", 1, 0)
cell_group_df$group2 <- ifelse(cell_group_df$celltype == "pos_scissor", 1, 0)

data1[["group1"]] <- cell_group_df$group1 # create cell names as metadata colum
data1[["group2"]] <- cell_group_df$group2 # create cell names as metadata colum
data1[["group3"]] <- ifelse(cell_group_df$group1 == 1, 2, 
												ifelse(cell_group_df$group2 == 1, 1, 0))

#table(data1[["group3"]])
#group3
#   0    1    2 
#1229  173  115 

###...07-3. remove scissor + cells
data_group <- subset(data1, idents = "pos_scissor", invert = TRUE)

###...07-4. plot of de genes 
my_vlnplot_for_neg_degene <- function(dt, gene) {
  VlnPlot(dt , features = gene , group.by = "group1",pt.size = 0, y.max = 8) +
  theme(legend.position = "top") +
  scale_fill_manual(values = c('grey','royalblue'), 
    labels = c("Others", "Scissor- cells")) +
  geom_signif(comparisons = list(c("0","1")),map_signif_level = TRUE, y_position =7.5) +
  theme(axis.title.x=element_blank()) + 
  scale_x_discrete(labels = c("",""))
}

my_umap_for_neg_degene <- function(dt, gene) {
  FeaturePlot(dt, features = gene)
}


# "SSR4" "RGS1" "HLA-DRB5" "APOE" "C1QB"  "C1QA"  "APOC1"  "JCHAIN"  "C1QC"  "DERL3"   

p1 <- my_vlnplot_for_neg_degene(data_group, "SSR4")
p2 <- my_vlnplot_for_neg_degene(data_group, "RGS1")
p3 <- my_vlnplot_for_neg_degene(data_group, "HLA-DRB5")
p4 <- my_vlnplot_for_neg_degene(data_group, "APOE")
p5 <- my_vlnplot_for_neg_degene(data_group, "C1QB")

p6 <- my_vlnplot_for_neg_degene(data_group, "C1QA")
p7 <- my_vlnplot_for_neg_degene(data_group, "APOC1")
p8 <- my_vlnplot_for_neg_degene(data_group, "JCHAIN")
p9 <- my_vlnplot_for_neg_degene(data_group, "C1QC")
p10 <- my_vlnplot_for_neg_degene(data_group, "DERL3")


# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

shared_legend <- extract_legend(p1)
grid.arrange(
p1 + theme(legend.position='hidden'),
p2 + theme(legend.position='hidden'),
p3 + theme(legend.position='hidden'),
p4 + theme(legend.position='hidden'),
p5 + theme(legend.position='hidden'),
p6 + theme(legend.position='hidden'),
p7 + theme(legend.position='hidden'),
p8 + theme(legend.position='hidden'),
p9 + theme(legend.position='hidden'),
p10 + theme(legend.position='hidden'), 
  top = shared_legend$grobs[[1]])


p1 <- my_umap_for_neg_degene(data1, "SSR4")
p2 <- my_umap_for_neg_degene(data1, "RGS1")
p3 <- my_umap_for_neg_degene(data1, "HLA-DRB5")
p4 <- my_umap_for_neg_degene(data1, "APOE")
p5 <- my_umap_for_neg_degene(data1, "C1QB")

p6 <- my_umap_for_neg_degene(data1, "C1QA")
p7 <- my_umap_for_neg_degene(data1, "APOC1")
p8 <- my_umap_for_neg_degene(data1, "JCHAIN")
p9 <- my_umap_for_neg_degene(data1, "C1QC")
p10 <- my_umap_for_neg_degene(data1, "DERL3")

shared_legend <- extract_legend(p4)
grid.arrange(
p1 + theme(legend.position='hidden'),
p2 + theme(legend.position='hidden'),
p3 + theme(legend.position='hidden'),
p4 + theme(legend.position='hidden'),
p5 + theme(legend.position='hidden'),
p6 + theme(legend.position='hidden'),
p7 + theme(legend.position='hidden'),
p8 + theme(legend.position='hidden'),
p9 + theme(legend.position='hidden'),
p10 + theme(legend.position='hidden'), 
  right = shared_legend$grobs[[1]])

    


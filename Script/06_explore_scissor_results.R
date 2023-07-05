#!/usr/bin/env Rscript
# 06_explore_scissor_results.R $INDIR 


### 01. Setting
source("Script/Functions/func_exploring.R")
args <- commandArgs(trailingOnly = TRUE)
library(Scissor)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
indir <- args[1]


### 02. Load data
load(paste0(indir, "infos1.RData"))
load(paste0(indir, "sc_dataset.RData"))



### 03. Extract Scissor cells
pos_cell_names <- infos1$Scissor_pos
neg_cell_names <- infos1$Scissor_neg



### 04. Check scissor cell distribution on umap plot 
data1 <- sc_dataset
Scissor_select <- rep(0, ncol(data1))
names(Scissor_select) <- colnames(data1)
Scissor_select[infos1$Scissor_pos] <- 1
Scissor_select[infos1$Scissor_neg] <- 2
data1 <- AddMetaData(data1, metadata = Scissor_select, col.name = "scissor")

# plot
p1 <- DimPlot(data1, reduction = 'umap', group.by = 'scissor', 
        cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
t <- length(Idents(sc_dataset))
n <- length(neg_cell_names)
p <- length(pos_cell_names)
labels_for_p <- c(paste0("Background cells (",t-p-n,")"),
                paste0("Scissor+ cell (",p,")"),
                paste0("Scissor- cell (",n,")"))

p1 <- p1 + scale_color_manual(values = c('grey','indianred1','royalblue'), 
					labels = labels_for_p)
p1 <- p1 + theme(legend.position = c(0.05, 0.2))


png(file = paste0(indir, "scissor_plot.png"))
p1
dev.off()




### 05. Plot for cluster of scissor cells
source("../Output/INTEGRATION/cell_type_ann_result/cluster_ann.R")

### ...05-1. make df for cluster of scissor cell ###
ndf <- make_cluster_n(neg_cell_names) # negative cells
pdf <- make_cluster_n(pos_cell_names) # positive cells

# merge negative and positive cells (long form)
ndf$type <- rep("Scissor- cells", nrow(ndf))
pdf$type <- rep("Scissor+ cells", nrow(pdf))
tdf <- rbind(ndf,pdf)
tdf$name <- factor(tdf$name, levels = paste0("cluster",  seq(0,14)))


### ...05-2. plot of cluster which consisted of scissor cells ###
# Define the colors for each cluster values for each type
cluster_colors <- 
c("cluster0" = "#f8766d", 
"cluster1" = "#E68613", 
"cluster2" = "#CD9600", 
"cluster3" = "#ABA300", 
"cluster4" = "#7CAE00", 
"cluster5" = "#0CB702", 
"cluster6" = "#00BE67",
"cluster7" = "#00C19A",
"cluster8" = "#00BFC4",
"cluster9" = "#00B8E7",
"cluster10" = "#00A9FF",
"cluster11" = "#C77CFF",
"cluster12" = "#ED68ED",
"cluster13" = "#FF61CC",
"cluster14" = "#FF68A1")

# Plot the bar plot
p2 <- ggplot(tdf, aes(name, n, fill = name)) +
  scale_y_break(c(60, 100)) +
  geom_bar(stat = 'identity') + facet_grid(. ~ type) +
  scale_fill_manual(values = cluster_colors) +
	geom_text(aes(label=ifelse(n == 0, "", n)), 
		position=position_dodge(width=0.9), vjust=-0.25) +
	theme_test()+ 
	xlab("") + ylab("The number of clusters")+
	theme(legend.position = "top",
		axis.text.x=element_text(angle=45, hjust=1),
		strip.text.x = element_text(size = 15),
		strip.background=element_rect(fill="white")) +
	guides(fill=guide_legend(title=""))


png(file = paste0(indir, "scissor_cluster.png"))
p2
dev.off()


### 06. Plot Which type of cell of scissor cell results ###
cell.name <- names(my_cluster_ann)
data1.new.idents <- case_when(
                      Idents(data1) %in% my_cluster_ann[[1]] ~  cell.name[1],
                      Idents(data1) %in% my_cluster_ann[[2]] ~  cell.name[2],
                      Idents(data1) %in% my_cluster_ann[[3]] ~  cell.name[3],
                      Idents(data1) %in% my_cluster_ann[[4]] ~  cell.name[4],
                      Idents(data1) %in% my_cluster_ann[[5]] ~  cell.name[5],
                      Idents(data1) %in% my_cluster_ann[[6]] ~  cell.name[6],
                      TRUE ~ NA_character_
                      )

###... 06-1. annotation clusters ###
# re-clustering for cluster 11
data2 <- subset(data1, idents = "11")
data2 <- RunPCA(data2, features= VariableFeatures(data2))
data2 <- FindNeighbors(object = data2, dims = 1:10)
data2 <- FindClusters(object = data2, resolution = 0.6)
data2 <- RunUMAP(object = data2, dims = 1:10)                    

bcell_names_in11 <- Idents(data2)[Idents(data2) == "0"] %>% names()                      
tcell_names_in11 <- Idents(data2)[Idents(data2) == "1"] %>% names()  

idxbcell <- names(Idents(data1)) %in% bcell_names_in11
idxtcell <- names(Idents(data1)) %in% tcell_names_in11
data1.new.idents[idxbcell] <- "B_cell"
data1.new.idents[idxtcell] <- "T_cell"
Idents(data1) <- data1.new.idents

data1[["CellName"]] <- colnames(data1) # create cell names as metadata colum
save(data1, file = paste0(indir, "annotated_data1.RData"))


###... 06-1. Subset Myeloid macrophage ###
# re-clustering for myeloid_macrophage
data3 <- subset(data1, idents = "Myeloid_macrophage")

selection.method = "vst"
resolution = 0.6
dims_Neighbors = 1:10
dims_TSNE = 1:10
dims_UMAP = 1:10
verbose = TRUE
DefaultAssay(data2) <- "RNA"
data3 <- FindVariableFeatures(object = data3 , selection.method = selection.method, verbose = verbose)
data3 <- ScaleData(object = data3, verbose = verbose)
data3 <- RunPCA(object = data3, features = VariableFeatures(data3), verbose = verbose)
data3 <- FindNeighbors(object = data3 , dims = dims_Neighbors, verbose = verbose)
data3 <- FindClusters( object = data3 , resolution = resolution, verbose = verbose)
data3 <- RunTSNE(object = data3 , dims = dims_TSNE)
data3 <- RunUMAP(object = data3 , dims = dims_UMAP, verbose = verbose)

### ...06-2. Check subpopulation of myeloid macrophage ###
DimPlot(data3, reduction = "umap") # subpopulation of myeloid marcrophage (5 clusters)

# https://www.cusabio.com/c-20938.html
# (source CD36)https://www.frontiersin.org/articles/10.3389/fimmu.2019.02035/full 
# (M1: CD86)  14.↵ P. J. Murray et al., Macrophage activation and polarization: nomenclature and experimental guidelines. Immunity 41, 14–20 (2014).CrossRefPubMedWeb of ScienceGoogle Scholar 
# M1 : CD40(https://www.nature.com/articles/s41598-021-87720-y)
 
M1marker = c("CD40", "CD86")
M2marker = c("CD163", "MS4A4A")
momarker = c("CD36")
VlnPlot(data3, features = momarker, pt.size = 0)
VlnPlot(data3, features = M1marker , pt.size = 0)
VlnPlot(data3, features = M2marker , pt.size = 0)
VlnPlot(data3, features = c(M1marker, M2marker, momarker) , pt.size = 0, ncol =2)
# cluster 0 : monocyte
# cluster 1,2,4 : M1/M2
# cluster 3 : M1

###...06-3. The umap visualization of myeloid macrophage subset ###
Scissor_select <- rep(0, ncol(data3))
names(Scissor_select) <- colnames(data3)
Scissor_select[infos1$Scissor_pos] <- 1
Scissor_select[infos1$Scissor_neg] <- 2
data2 <- AddMetaData(data3, metadata = Scissor_select, col.name = "scissor")

# The umap visualization of myeloid marcrophage subset
t <- Scissor_select[Scissor_select == 0] %>% length()
n <- Scissor_select[Scissor_select == 2] %>% length()
p <- Scissor_select[Scissor_select == 1] %>% length()
labels_for_p <- c(paste0("Background cells (",t,")"),
                paste0("Scissor+ cell (",p,")"),
                paste0("Scissor- cell (",n,")"))
DimPlot(data3, reduction = 'umap', group.by = 'scissor', 
        cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1)) +
	scale_color_manual(values = c('grey','indianred1','royalblue'), 
	 labels =labels_for_p) + 
	theme(legend.position = c(0.6, 0.2))

###...06-4. make df for barplot of cell type ###
cellndf <- make_celltype_n(neg_cell_names) # negative cells
cellpdf <- make_celltype_n(pos_cell_names) # positive cells

# merge negative and positive cells (long form)
cellndf$type <- rep("Scissor- cells", nrow(cellndf))
cellpdf$type <- rep("Scissor+ cells", nrow(cellpdf))
celltdf <- rbind(cellndf,cellpdf)

celltdf[celltdf$type == "Scissor+ cells" & celltdf$name == "Myeloid macrophage",]$name <- "Monocyte"
celltdf[celltdf$type == "Scissor- cells" & celltdf$name == "Myeloid macrophage",]$name <- "M1/M2 macrophage"
celltdf$name <- factor(celltdf$name, levels = c("M1/M2 macrophage","Monocyte", "B cell", "Fibroblast",
		"T cell", "other"))

###... 06-5. plot for cell type ###
celltype_colors <- 
c("M1/M2 macrophage" = "#fd4527", 
"Monocyte" = "#fda927",
"B cell" = "#439e27", 
"Fibroblast" = "#439ea5", 
"T cell" = "#8a5098", 
"other" = "gray30")

p3 <- ggplot(celltdf, aes(name, No, fill = name)) +
  geom_bar(stat = 'identity') + facet_grid(. ~ type) +
  scale_fill_manual(values = celltype_colors ) +
	geom_text(aes(label=ifelse(No == 0, "", No)), 
		position=position_dodge(width=0.9), vjust=-0.25) +
	theme_test()+ 
	xlab("") + ylab("The number of cell type")+
	theme(legend.position = "top",
		axis.text.x=element_text(angle=45, hjust=1),
		strip.text.x = element_text(size = 15),
		strip.background=element_rect(fill="white")) +
	guides(fill=guide_legend(title=""))


png(file = paste0(indir, "scissor_celltype.png"))
p3
dev.off()







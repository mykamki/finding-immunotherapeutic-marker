#!/usr/bin/env Rscript
# 06_explore_scissor_results.R $INDIR 


### 01. Setting
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




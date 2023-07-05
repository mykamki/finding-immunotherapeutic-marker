#!/usr/bin/env Rscript
# 06_explore_scissor_results.R $INDIR 

#### 01.
args <- commandArgs(trailingOnly = TRUE)
library(Scissor)
library(ggplot2)
library(dplyr)
library(RColorBrewer)


#### 01. Load data ####
location <- args[1]
load(paste0(location, "infos1.RData"))
load(paste0(location, "preprocessed_data1.RData"))



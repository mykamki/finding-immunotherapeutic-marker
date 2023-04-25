#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(cluster, quietly = TRUE)
library(dplyr)
library(boot)
library(ggplot2)
indir <- args[1]



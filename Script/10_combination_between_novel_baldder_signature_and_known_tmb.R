### 01. Setting
library(survival)
library(dplyr)
library(survminer)
library(cowplot)
library(grid)
library(gridExtra)
library(ggplot2)
library(reshape2)
options(scipen=999)
indir <- args[1]
outdir <- args[2]
source("Script/Functions/func_clinical_relevance.R")



### 02. Load data
load(paste0(indir,"clinical_gse176307.RData"))
load(paste0(indir,"clinical_imvigor210core.RData"))
load(paste0(indir,"clinical_ucgenome.RData"))


### 03. Survival relevance and ICB response of combination between novel bladder signature and known TMB in GSE176307
clinical_gse176307_2 <- make_combination_group_dataset(clinical_gse176307)
clinical_gse176307_2 <- make_binary_response_dataset(clinical_gse176307_2)
p1a <- my_combination_survival_plot(clinical_gse176307_2, "GSE176307")
p1b <- bar_plot_by_combination(clinical_gse176307_2)


### 04. Survival relevance and ICB response of ncombination between novel bladder signature and known TMB in IMvigor210
clinical_imvigor210core_2 <- make_combination_group_dataset(clinical_imvigor210core)
clinical_imvigor210core_2 <- make_binary_response_dataset(clinical_imvigor210core_2)
p2a <- my_combination_survival_plot(clinical_imvigor210core_2, "IMvigor210")
p2b <- bar_plot_by_combination(clinical_imvigor210core_2)


### 05. Survival relevance and ICB response of combination between novel bladder signature and known TMB in UC-GENOME
clinical_ucgenome_2 <- make_combination_group_dataset(clinical_ucgenome)
clinical_ucgenome_2 <- make_binary_response_dataset(clinical_ucgenome_2)
p3a <- my_combination_survival_plot(clinical_ucgenome_2, "UC-GENOME")
p3b <- bar_plot_by_combination(clinical_ucgenome_2)


### 06. Plotting
pA <- grid.arrange(
grid.arrange(p1a$plot + theme(legend.position='hidden'),p1a$table, layout_matrix = rbind(c(1), c(1), c(2))),
p1b + theme(legend.position='hidden'),
grid.arrange(p2a$plot + theme(legend.position='hidden'),p2a$table, layout_matrix = rbind(c(1), c(1), c(2))),
p2b + theme(legend.position='hidden'),
grid.arrange(p3a$plot + theme(legend.position='hidden'),p3a$table, layout_matrix = rbind(c(1), c(1), c(2))),
p3b + theme(legend.position='hidden'),
ncol = 2)

bar_legend <- extract_legend(p1b)
sur1_legend <- extract_legend(p1a$plot)
shared_legend <- grid.arrange(sur1_legend, bar_legend, nrow =1)

pB <- grid.arrange(
grid.arrange(p1a$plot + theme(legend.position='hidden'),p1a$table, layout_matrix = rbind(c(1), c(1), c(2))),
grid.arrange(p2a$plot + theme(legend.position='hidden'),p2a$table, layout_matrix = rbind(c(1), c(1), c(2))),
grid.arrange(p3a$plot + theme(legend.position='hidden'),p3a$table, layout_matrix = rbind(c(1), c(1), c(2))),
ncol = 1)


### 07. Multiple plots
source("Script/Functions/hh_vs_other.R")
source("Script/Functions/other_vs_ll.R")
source("Script/Functions/make_all_survival_plot.R")


com4 <- grid.arrange(pA, shared_legend, heights = c(8, 1))
com4 <- grid.arrange(pB, sur1_legend, heights = c(8, 1))
comper <- make_all_survival_plot(clinical_gse176307, clinical_imvigor210core, clinical_ucgenome)

phh <- hh_vs_other(clinical_gse176307, clinical_imvigor210core, clinical_ucgenome)
hh <- my_anno_plot("Both High vs. Any Low")

pll <- other_vs_ll(clinical_gse176307, clinical_imvigor210core, clinical_ucgenome)
ll <- my_anno_plot("Any High vs. Both Low")


# hh vs. other / other vs. ll
png(paste0(outdir, "Two_group_combination.png"),  width = 800, height = 900)
grid.arrange(ggarrange(hh,ll,nrow =1), ggarrange(phh, pll, nrow=1), heights = c(0.5, 8))
dev.off()

# all survival plots
png(paste0(outdir, "all_sur_plot.png"),  width = 900, height = 900)
grid.arrange(com4, comper, widths =c (1.2,3))
dev.off()

#p <- grid.arrange(qa, qb, qc, shared_legend, heights = c(2,2,2,0.5))
#png(paste0(outdir, "supfig6.png"),  width = 600, height = 750)
#as_ggplot(p) +  draw_plot_label(label = c("a", "b", "c"), size = 15,
#                  x = c(0, 0, 0), y = c(1, 0.7, 0.4))   # transform to a ggplot
#dev.off()

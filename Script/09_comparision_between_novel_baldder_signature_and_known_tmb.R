
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


### 03. Survival relevance and ICB response of novel bladder signature in GSE176307
p1a <- my_geneset_survival_plot(clinical_gse176307, "GSE176307")
clinical_gse176307 <- make_binary_response_dataset(clinical_gse176307)
p1b <- bar_plot_by_gene(clinical_gse176307)


### 04. Survival relevance and ICB response of novel bladder signature in IMvigor210
p2a <- my_geneset_survival_plot(clinical_imvigor210core, "IMvigor210")
clinical_imvigor210core <- make_binary_response_dataset(clinical_imvigor210core)
p2b <- bar_plot_by_gene(clinical_imvigor210core)


### 05. Survival relevance and ICB response of novel bladder signature in UC-GENOME
p3a <- my_geneset_survival_plot(clinical_ucgenome, "UC-GENOME")
clinical_ucgenome <- make_binary_response_dataset(clinical_ucgenome)
p3b <- bar_plot_by_gene(clinical_ucgenome)


### 06. Survival relevance and ICB response of known TMB in GSE176307
clinical_gse176307_2 <- make_tmb_group_dataset(clinical_gse176307)
q1a <- my_tmb_survival_plot(clinical_gse176307_2, "GSE176307")
q1b <- bar_plot_by_tmb(clinical_gse176307_2)


### 07. Survival relevance and ICB response of known TMB in IMvigor210
clinical_imvigor210core_2 <- make_tmb_group_dataset(clinical_imvigor210core)
q2a <- my_tmb_survival_plot(clinical_imvigor210core_2, "IMvigor210")
q2b <- bar_plot_by_tmb(clinical_imvigor210core_2)


### 08. Survival relevance and ICB response of known TMB in UC-GENOME
clinical_ucgenome_2 <- make_tmb_group_dataset(clinical_ucgenome)
q3a <- my_tmb_survival_plot(clinical_ucgenome_2, "UC-GENOME")
q3b <- bar_plot_by_tmb(clinical_ucgenome_2)


### 09. Plotting
pA <- grid.arrange(
grid.arrange(p1a$plot + theme(legend.position='hidden'),p1a$table, layout_matrix = rbind(c(1), c(1), c(2))),
p1b + theme(legend.position='hidden'),
grid.arrange(p2a$plot + theme(legend.position='hidden'),p2a$table, layout_matrix = rbind(c(1), c(1), c(2))),
p2b + theme(legend.position='hidden'),
grid.arrange(p3a$plot + theme(legend.position='hidden'),p3a$table, layout_matrix = rbind(c(1), c(1), c(2))),
p3b + theme(legend.position='hidden'),
ncol = 2)

pB <- grid.arrange(
grid.arrange(q1a$plot + theme(legend.position='hidden'),q1a$table, layout_matrix = rbind(c(1), c(1), c(2))),
q1b + theme(legend.position='hidden'),
grid.arrange(q2a$plot + theme(legend.position='hidden'),q2a$table, layout_matrix = rbind(c(1), c(1), c(2))),
q2b + theme(legend.position='hidden'),
grid.arrange(q3a$plot + theme(legend.position='hidden'),q3a$table, layout_matrix = rbind(c(1), c(1), c(2))),
q3b + theme(legend.position='hidden'),
ncol = 2)

bar_legend <- extract_legend(q3b)
sur1_legend <- extract_legend(p1a$plot)
sur2_legend <- extract_legend(q1a$plot)
shared_legend <- grid.arrange(sur1_legend, sur2_legend, bar_legend, nrow =1)

anno1 <- my_anno_plot("Novel Gene Signature")
anno2 <- my_anno_plot("Known TMB")

png(paste0(outdir, "Novel_signature_and_tmb.png"),  width = 800, height = 800)
grid.arrange(ggarrange(anno1, anno2, nrow =1), ggarrange(pA,pB),shared_legend, heights = c(0.5,8,0.5))
dev.off()


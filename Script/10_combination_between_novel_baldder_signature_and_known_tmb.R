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

my_combination_survival_plot <- function(dataset, datasettitle) {
	fit1 <- survfit(Surv(time = time , event = status )~ Combination_Signature_TMB , data = dataset)
	ggsurvplot(fit1, data = dataset, pval = T, pval.size =6, risk.table= T, 
                   title= datasettitle, legend = "right",
		   #palette = c("#00BFC4", "#F8766D"),
                   ggtheme=custom_theme(), tables.height = 0.1,
                   tables.theme = theme_cleantable(), fontsize = 3, 
                   risk.table.y.text.col = T, # colour risk table text annotations.
                   risk.table.y.text = F, # show bars instead of names in text annotations
                   #legend.labs = expression("Signature[High]", "Signature[Low]")  # Change legend labels
                )
}

my_extreme_combination_survival_plot <- function(dataset, datasettitle) {
	dataset <- dataset %>% filter(Combination_Signature_TMB %in% c("HL", "LL"))
	fit1 <- survfit(Surv(time = time , event = status )~ Combination_Signature_TMB , data = dataset)
	ggsurvplot(fit1, data = dataset, pval = T, pval.size =6, risk.table= T, 
                   title= datasettitle, legend = "right",
                   ggtheme=custom_theme(), tables.height = 0.1,
                   tables.theme = theme_cleantable(), fontsize = 3, 
                   risk.table.y.text.col = T, # colour risk table text annotations.
                   risk.table.y.text = F, # show bars instead of names in text annotations
                   #legend.labs = expression("Signature[High]", "Signature[Low]")  # Change legend labels
                )
}

make_combination_group_dataset <- function(dataset) {
	dataset <- dataset %>% filter(!is.na(tmb)) 
	tmb_med <- dataset %>% select(tmb) %>% as.matrix() %>% as.vector() %>% median()
	dataset <- dataset %>% mutate(Known_TMB = ifelse(tmb >= tmb_med, "High", 
							 ifelse(tmb < tmb_med, "Low", NA)))
	dataset <- dataset %>% mutate(
		Combination_Signature_TMB = case_when(
			Novel_Signature == "High" & Known_TMB == "High" ~ "HH",
			Novel_Signature == "High" & Known_TMB == "Low" ~ "HL",
			Novel_Signature == "Low" & Known_TMB == "High" ~ "LH",
			Novel_Signature == "Low" & Known_TMB == "Low" ~ "LL",
			TRUE ~ NA_character_
			)
		)
	#dataset$Combination_Signature_TMB <- factor(dataset$Combination_Signature_TMB, levels = c("Both High", "Any Low"))
	return(dataset)
}

make_combination_group_dataset <- function(dataset) {
	dataset <- dataset %>% filter(!is.na(tmb)) 
	tmb_med <- dataset %>% select(tmb) %>% as.matrix() %>% as.vector() %>% median()
	dataset <- dataset %>% mutate(Known_TMB = ifelse(tmb >= tmb_med, "High", 
							 ifelse(tmb < tmb_med, "Low", NA)))
	dataset <- dataset %>% mutate(
		Combination_Signature_TMB = case_when(
			Novel_Signature == "High" & Known_TMB == "High" ~ "HH",
			Novel_Signature == "High" & Known_TMB == "Low" ~ "Any",
			Novel_Signature == "Low" & Known_TMB == "High" ~ "Any",
			Novel_Signature == "Low" & Known_TMB == "Low" ~ "Any",
			TRUE ~ NA_character_
			)
		)
	dataset$Combination_Signature_TMB <- factor(dataset$Combination_Signature_TMB, levels = c("HH", "Any"))
	return(dataset)
}

make_combination_group_dataset <- function(dataset) {
	dataset <- dataset %>% filter(!is.na(tmb)) 
	tmb_med <- dataset %>% select(tmb) %>% as.matrix() %>% as.vector() %>% median()
	dataset <- dataset %>% mutate(Known_TMB = ifelse(tmb >= tmb_med, "High", 
							 ifelse(tmb < tmb_med, "Low", NA)))
	dataset <- dataset %>% mutate(
		Combination_Signature_TMB = case_when(
			Novel_Signature == "High" & Known_TMB == "High" ~ "Any",
			Novel_Signature == "High" & Known_TMB == "Low" ~ "Any",
			Novel_Signature == "Low" & Known_TMB == "High" ~ "Any",
			Novel_Signature == "Low" & Known_TMB == "Low" ~ "LL",
			TRUE ~ NA_character_
			)
		)
	#dataset$Combination_Signature_TMB <- factor(dataset$Combination_Signature_TMB, levels = c("Both High", "Any Low"))
	return(dataset)
}

bar_plot_by_combination <- function(dataset) {
	data <- dataset %>% 
                filter(!is.na(response.binary)) %>% 
                group_by(Combination_Signature_TMB, response.binary) %>% 
                summarize(n = n())
        data$response.binary <- factor(data$response.binary, levels = c("R", "NR"))
	
	hh <- sum(data[data$Combination_Signature_TMB == "HH",]$n)
	hl <- sum(data[data$Combination_Signature_TMB == "HL",]$n)
	lh <- sum(data[data$Combination_Signature_TMB == "LH",]$n)
	ll <- sum(data[data$Combination_Signature_TMB == "LL",]$n)

	data$tt <- c(rep(hh,nrow(data[data$Combination_Signature_TMB == "HH",])), 
		     rep(hl,nrow(data[data$Combination_Signature_TMB == "HL",])), 
		     rep(lh,nrow(data[data$Combination_Signature_TMB == "LH",])), 
		     rep(ll,nrow(data[data$Combination_Signature_TMB == "LL",])))
	data <- data %>% mutate(percent = round(n/tt*100,1)) 

	dt <- dcast(data, Combination_Signature_TMB ~ response.binary, value.var="n")
	change_na2zero <- function(x) {ifelse(is.na(x), 0, as.double(x))}
	apply(dt[,2:3], 2, change_na2zero) %>% fisher.test(simulate.p.value=TRUE,B=1e7) -> dtt

        # make plot
        p <- ggplot(data, aes(fill=response.binary, y=percent, x=Combination_Signature_TMB)) +
                geom_col() + 
                geom_signif(comparisons = list(c("HH", "LL")), 
                            annotations = paste0("pvalue=", round(dtt$p.value,2)), y_position = 102) +
                scale_x_discrete(labels=c(paste0("HH\n(n=", hh,")"),
					  paste0("HL\n(n=", hl,")"), 
					  paste0("LH\n(n=", lh,")"), 
					  paste0("LL\n(n=", ll,")"))) +
                scale_fill_manual(values = c("orange", "royalblue3")) +
                xlab("") + ylab("") + 
                theme(panel.background = element_blank(),
                     axis.ticks = element_blank(),
                     axis.text.y = element_blank(),legend.position="top",
                     plot.margin = margin(0, 1, 0, 1, "cm")) +
                geom_text(aes(label=paste0(n," (",percent,"%)")), position = position_stack(vjust = 0.5)) +
                guides(fill=guide_legend(title="Responder")) + 
                theme(axis.text.x = element_text(size = 12))
	
	return(p)
}

bar_plot_by_combination <- function(dataset) {
	#dataset <- dataset %>% filter(Combination_Signature_TMB %in% c("HH", "LL"))
	data <- dataset %>% 
                filter(!is.na(response.binary)) %>% 
                group_by(Combination_Signature_TMB, response.binary) %>% 
                summarize(n = n())
        data$response.binary <- factor(data$response.binary, levels = c("R", "NR"))
	
	hh <- sum(data[data$Combination_Signature_TMB == "Any",]$n)
	ll <- sum(data[data$Combination_Signature_TMB == "LL",]$n)

	data$tt <- c(rep(hh,nrow(data[data$Combination_Signature_TMB == "Any",])), 
		     rep(ll,nrow(data[data$Combination_Signature_TMB == "LL",])))
	data <- data %>% mutate(percent = round(n/tt*100,1)) 

	dt <- dcast(data, Combination_Signature_TMB ~ response.binary, value.var="n")
	change_na2zero <- function(x) {ifelse(is.na(x), 0, as.double(x))}
	apply(dt[,2:3], 2, change_na2zero) %>% fisher.test(simulate.p.value=TRUE,B=1e7) -> dtt

        # make plot
        p <- ggplot(data, aes(fill=response.binary, y=percent, x=Combination_Signature_TMB)) +
                geom_col() + 
                geom_signif(comparisons = list(c("Any", "LL")), 
                            annotations = paste0("pvalue=", round(dtt$p.value,2)), y_position = 102) +
                scale_x_discrete(labels=c(paste0("Any\n(n=", hh,")"), 
					  paste0("LL\n(n=", ll,")"))) +
                scale_fill_manual(values = c("orange", "royalblue3")) +
                xlab("") + ylab("") + 
                theme(panel.background = element_blank(),
                     axis.ticks = element_blank(),
                     axis.text.y = element_blank(),legend.position="top",
                     plot.margin = margin(0, 1, 0, 1, "cm")) +
                geom_text(aes(label=paste0(n," (",percent,"%)")), position = position_stack(vjust = 0.5)) +
                guides(fill=guide_legend(title="Responder")) + 
                theme(axis.text.x = element_text(size = 12))
	
	return(p)
}

make_combination_group_dataset <- function(dataset) {
        dataset <- dataset %>% filter(!is.na(tmb)) 
        tmb_med <- dataset %>% select(tmb) %>% as.matrix() %>% as.vector() %>% median()
        dataset <- dataset %>% mutate(Known_TMB = ifelse(tmb >= tmb_med, "High", 
                                                         ifelse(tmb < tmb_med, "Low", NA)))
        dataset <- dataset %>% mutate(
                Combination_Signature_TMB = case_when(
                        Novel_Signature == "High" & Known_TMB == "High" ~ "HH",
                        Novel_Signature == "High" & Known_TMB == "Low" ~ "HL",
                        Novel_Signature == "Low" & Known_TMB == "High" ~ "LH",
                        Novel_Signature == "Low" & Known_TMB == "Low" ~ "LL",
                        TRUE ~ NA_character_
                        )
                )
        #dataset$Combination_Signature_TMB <- factor(dataset$Combination_Signature_TMB, levels = c("Both High", "Any Low"))
        return(dataset)
}

### 02. Load data
load(paste0(indir,"clinical_gse176307.RData"))
load(paste0(indir,"clinical_imvigor210core.RData"))
load(paste0(indir,"clinical_ucgenome.RData"))


### 03. Survival relevance and ICB response of combination between novel bladder signature and known TMB in GSE176307
clinical_gse176307 <- make_combination_group_dataset(clinical_gse176307)
clinical_gse176307 <- make_binary_response_dataset(clinical_gse176307)
p1a <- my_combination_survival_plot(clinical_gse176307, "GSE176307")
p1b <- bar_plot_by_combination(clinical_gse176307)


### 04. Survival relevance and ICB response of ncombination between novel bladder signature and known TMB in IMvigor210
clinical_imvigor210core <- make_combination_group_dataset(clinical_imvigor210core)
clinical_imvigor210core <- make_binary_response_dataset(clinical_imvigor210core)
p2a <- my_combination_survival_plot(clinical_imvigor210core, "IMvigor210")
p2b <- bar_plot_by_combination(clinical_imvigor210core)


### 05. Survival relevance and ICB response of combination between novel bladder signature and known TMB in UC-GENOME
clinical_ucgenome <- make_combination_group_dataset(clinical_ucgenome)
clinical_ucgenome <- make_binary_response_dataset(clinical_ucgenome)
p3a <- my_combination_survival_plot(clinical_ucgenome, "UC-GENOME")
p3b <- bar_plot_by_combination(clinical_ucgenome)



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
			
png(paste0(outdir, "fig3.png"),  width = 500, height = 800)
grid.arrange(pA, shared_legend, heights = c(8, 1))
dev.off()


png(paste0(outdir, "test.png"),  width = 400, height = 800)
grid.arrange(pA, shared_legend, heights = c(8, 1))
dev.off()

my_HH_survival_plot <- function(dataset, datasettitle) {
	dataset <- dataset %>% filter(Combination_Signature_TMB %in% c("HH", "LL"))
	fit1 <- survfit(Surv(time = time , event = status )~ Combination_Signature_TMB , data = dataset)
	ggsurvplot(fit1, data = dataset, pval = T, pval.size =6, risk.table= T, 
                   title= datasettitle, legend = "right",
		   palette = c("#F8766D", "#C77CFF"),
                   ggtheme=custom_theme(), tables.height = 0.1,
                   tables.theme = theme_cleantable(), fontsize = 3, 
                   risk.table.y.text.col = T, # colour risk table text annotations.
                   risk.table.y.text = F, # show bars instead of names in text annotations
                   #legend.labs = expression("Signature[High]", "Signature[Low]")  # Change legend labels
                )
}

my_LH_survival_plot <- function(dataset, datasettitle) {
	dataset <- dataset %>% filter(Combination_Signature_TMB %in% c("LH", "LL"))
	fit1 <- survfit(Surv(time = time , event = status )~ Combination_Signature_TMB , data = dataset)
	ggsurvplot(fit1, data = dataset, pval = T, pval.size =6, risk.table= T, 
                   title= datasettitle, legend = "right",
		   palette = c("#00BFC4", "#C77CFF"),
                   ggtheme=custom_theme(), tables.height = 0.1,
                   tables.theme = theme_cleantable(), fontsize = 3, 
                   risk.table.y.text.col = T, # colour risk table text annotations.
                   risk.table.y.text = F, # show bars instead of names in text annotations
                   #legend.labs = expression("Signature[High]", "Signature[Low]")  # Change legend labels
                )
}

my_HL_survival_plot <- function(dataset, datasettitle) {
	dataset <- dataset %>% filter(Combination_Signature_TMB %in% c("HL", "LL"))
	fit1 <- survfit(Surv(time = time , event = status )~ Combination_Signature_TMB , data = dataset)
	ggsurvplot(fit1, data = dataset, pval = T, pval.size =6, risk.table= T, 
                   title= datasettitle, legend = "right",
		   palette = c("#7CAE00", "#C77CFF"),
                   ggtheme=custom_theme(), tables.height = 0.1,
                   tables.theme = theme_cleantable(), fontsize = 3, 
                   risk.table.y.text.col = T, # colour risk table text annotations.
                   risk.table.y.text = F, # show bars instead of names in text annotations
                   #legend.labs = expression("Signature[High]", "Signature[Low]")  # Change legend labels
                )
}

qa1 <- my_HH_survival_plot(clinical_gse176307, "")
qa2 <- my_LH_survival_plot(clinical_gse176307, "")
qa3 <- my_HL_survival_plot(clinical_gse176307, "")

qb1 <- my_HH_survival_plot(clinical_imvigor210core, "")
qb2 <- my_LH_survival_plot(clinical_imvigor210core, "")
qb3 <- my_HL_survival_plot(clinical_imvigor210core, "")

qc1 <- my_HH_survival_plot(clinical_ucgenome, "")
qc2 <- my_LH_survival_plot(clinical_ucgenome, "")
qc3 <- my_HL_survival_plot(clinical_ucgenome, "")


qa <- grid.arrange(
grid.arrange(qa1$plot + theme(legend.position='hidden'),qa1$table, layout_matrix = rbind(c(1), c(1), c(2))),
grid.arrange(qa2$plot + theme(legend.position='hidden'),qa2$table, layout_matrix = rbind(c(1), c(1), c(2))),
grid.arrange(qa3$plot + theme(legend.position='hidden'),qa3$table, layout_matrix = rbind(c(1), c(1), c(2))),
nrow = 1)

qb <- grid.arrange(
grid.arrange(qb1$plot + theme(legend.position='hidden'),qb1$table, layout_matrix = rbind(c(1), c(1), c(2))),
grid.arrange(qb2$plot + theme(legend.position='hidden'),qb2$table, layout_matrix = rbind(c(1), c(1), c(2))),
grid.arrange(qb3$plot + theme(legend.position='hidden'),qb3$table, layout_matrix = rbind(c(1), c(1), c(2))),
nrow = 1)

qc <- grid.arrange(
grid.arrange(qc1$plot + theme(legend.position='hidden'),qc1$table, layout_matrix = rbind(c(1), c(1), c(2))),
grid.arrange(qc2$plot + theme(legend.position='hidden'),qc2$table, layout_matrix = rbind(c(1), c(1), c(2))),
grid.arrange(qc3$plot + theme(legend.position='hidden'),qc3$table, layout_matrix = rbind(c(1), c(1), c(2))),
nrow = 1)

sur1_legend <- extract_legend(qa1$plot)
sur2_legend <- extract_legend(qa2$plot)
sur3_legend <- extract_legend(qa3$plot)

shared_legend <- grid.arrange(sur1_legend, sur2_legend, sur3_legend, nrow =1)


p <- grid.arrange(qa, qb, qc, shared_legend, heights = c(2,2,2,0.5))

png(paste0(outdir, "supfig6.png"),  width = 600, height = 750)
as_ggplot(p) +  draw_plot_label(label = c("a", "b", "c"), size = 15,
                  x = c(0, 0, 0), y = c(1, 0.7, 0.4))   # transform to a ggplot
dev.off()

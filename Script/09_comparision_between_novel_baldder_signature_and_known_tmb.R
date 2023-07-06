
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


custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}

my_geneset_survival_plot <- function(dataset, datasettitle) {
	fit1 <- survfit(Surv(time = time , event = status )~ Novel_Signature , data = dataset)
	ggsurvplot(fit1, data = dataset, pval = T, pval.size =6, risk.table= T, 
                   title= datasettitle, 
                   ggtheme=custom_theme(), tables.height = 0.1,
                   tables.theme = theme_cleantable(), fontsize = 3, 
                   risk.table.y.text.col = T, # colour risk table text annotations.
                   risk.table.y.text = F, # show bars instead of names in text annotations
                   legend.labs = expression("Signature[High]", "Signature[Low]")  # Change legend labels
                )
}

make_tmb_group_dataset <- function(dataset) {
	dataset <- dataset %>% filter(!is.na(tmb)) 
	tmb_med <- dataset %>% select(tmb) %>% as.matrix() %>% as.vector() %>% median()
	dataset <- dataset %>% mutate(Known_TMB = ifelse(tmb >= tmb_med, "High", "Low"))
	return(dataset)
}

my_tmb_survival_plot <- function(dataset, datasettitle) {
	fit1 <- survfit(Surv(time = time , event = status )~ Known_TMB , data = dataset)
	ggsurvplot(fit1, data = dataset, pval = T, pval.size =6, risk.table= T, 
                   title= datasettitle, 
                   ggtheme=custom_theme(), tables.height = 0.1,
                   tables.theme = theme_cleantable(), fontsize = 3, 
                   risk.table.y.text.col = T, # colour risk table text annotations.
                   risk.table.y.text = F, # show bars instead of names in text annotations
                   legend.labs = expression("Signature[High]", "Signature[Low]")  # Change legend labels
                )
}

make_binary_response_dataset <- function(dataset) {
  dataset$response.binary <- ifelse(dataset$response %in% c("CR", "PR"), "R", 
                              ifelse(dataset$response %in% c("SD", "PD"), "NR", NA))
  return(dataset)
  }

bar_plot_by_gene <- function(dataset) {
	data <- dataset %>% 
                filter(!is.na(response.binary)) %>% 
                group_by(Novel_Signature, response.binary) %>% 
                summarize(n = n())
        data$response.binary <- factor(data$response.binary, levels = c("R", "NR"))
	
	ht <- sum(data[data$Novel_Signature == "High",]$n)
        lt <- sum(data[data$Novel_Signature == "Low",]$n)
        data$tt <- c(rep(ht, nrow(data[data$Novel_Signature == "High",])), rep(lt,nrow(data[data$Novel_Signature == "Low",])))
        data <- data %>% mutate(percent = round(n/tt*100,0)) 

	dt <- dcast(data, Novel_Signature ~ response.binary, value.var="n")
        change_na2zero <- function(x) {ifelse(is.na(x), 0, as.double(x))}
        apply(dt[,2:3], 2, change_na2zero) %>% fisher.test() -> dtt

        # make plot
        p <- ggplot(data, aes(fill=response.binary, y=percent, x=Novel_Signature)) +
                geom_col() + 
                geom_signif(comparisons = list(c("High", "Low")), 
                            annotations = paste0("pvalue=", round(dtt$p.value,2)), y_position = 102) +
                scale_x_discrete(labels=c(paste0("High\n(n=", ht,")"), paste0("Low\n(n=", lt,")"))) +
                scale_fill_manual(values = c("orange", "royalblue3")) +
                xlab("") + ylab("") + 
                theme(panel.background = element_blank(),
                     axis.ticks = element_blank(),
                     axis.text.y = element_blank(),
                     plot.margin = margin(0, 1, 0, 1, "cm")) +
                geom_text(aes(label=paste0(n," (",percent,"%)")), position = position_stack(vjust = 0.5)) +
                guides(fill=guide_legend(title="Responder")) + 
                theme(axis.text.x = element_text(size = 12))
	
	return(p)
}

bar_plot_by_tmb <- function(dataset) {
	data <- dataset %>% 
                filter(!is.na(response.binary)) %>% 
                group_by(Known_TMB, response.binary) %>% 
                summarize(n = n())
        data$response.binary <- factor(data$response.binary, levels = c("R", "NR"))
	
	ht <- sum(data[data$Known_TMB == "High",]$n)
        lt <- sum(data[data$Known_TMB == "Low",]$n)
        data$tt <- c(rep(ht, nrow(data[data$Known_TMB == "High",])), rep(lt,nrow(data[data$Known_TMB == "Low",])))
        data <- data %>% mutate(percent = round(n/tt*100,0)) 

	dt <- dcast(data, Known_TMB ~ response.binary, value.var="n")
        change_na2zero <- function(x) {ifelse(is.na(x), 0, as.double(x))}
        apply(dt[,2:3], 2, change_na2zero) %>% fisher.test() -> dtt

        # make plot
        p <- ggplot(data, aes(fill=response.binary, y=percent, x=Known_TMB)) +
                geom_col() + 
                geom_signif(comparisons = list(c("High", "Low")), 
                            annotations = paste0("pvalue=", round(dtt$p.value,2)), y_position = 102) +
                scale_x_discrete(labels=c(paste0("High\n(n=", ht,")"), paste0("Low\n(n=", lt,")"))) +
                scale_fill_manual(values = c("orange", "royalblue3")) +
                xlab("") + ylab("") + 
                theme(panel.background = element_blank(),
                     axis.ticks = element_blank(),
                     axis.text.y = element_blank(),
                     plot.margin = margin(0, 1, 0, 1, "cm")) +
                geom_text(aes(label=paste0(n," (",percent,"%)")), position = position_stack(vjust = 0.5)) +
                guides(fill=guide_legend(title="Responder")) + 
                theme(axis.text.x = element_text(size = 12))
	
	return(p)
}


### 02. Load data
load(paste0(indir,"clinical_gse176307.RData"))
load(paste0(indir,"clinical_imvigor210core.RData"))
load(paste0(indir,"clinical_ucgenome.RData"))


pA1 <- my_geneset_survival_plot(clinical_gse176307, "GSE176307")
pA2 <- my_geneset_survival_plot(clinical_imvigor210core, "IMvigor210")
pA3 <- my_geneset_survival_plot(clinical_ucgenome, "UC-GENOME")

png("surv.png")
pA <- grid.arrange(
grid.arrange(pA1$plot + theme(legend.position='hidden'),pA1$table, layout_matrix = rbind(c(1), c(1), c(2))),
grid.arrange(pA2$plot + theme(legend.position='hidden'),pA2$table, layout_matrix = rbind(c(1), c(1), c(2))),
grid.arrange(pA3$plot + theme(legend.position='hidden'),pA3$table, layout_matrix = rbind(c(1), c(1), c(2))),
ncol = 3)
dev.off()

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
q1a <- my_tmb_survival_plot(clinical_gse176307, "GSE176307")
clinical_gse176307 <- make_tmb_group_dataset(clinical_gse176307)
q1b <- bar_plot_by_tmb(clinical_gse176307)


### 07. Survival relevance and ICB response of known TMB in IMvigor210
q2a <- my_tmb_survival_plot(clinical_imvigor210core, "IMvigor210")
clinical_imvigor210core <- make_tmb_group_dataset(clinical_imvigor210core)
q2b <- bar_plot_by_tmb(clinical_imvigor210core)


### 08. Survival relevance and ICB response of known TMB in UC-GENOME
q3a <- my_tmb_survival_plot(clinical_ucgenome, "UC-GENOME")
clinical_ucgenome <- make_tmb_group_dataset(clinical_ucgenome)
q3b <- bar_plot_by_tmb(clinical_ucgenome)


### 09. 
grid.arrange(
grid.arrange(p1a$plot + theme(legend.position='hidden'),p1a$table, layout_matrix = rbind(c(1), c(1), c(2))),
p1b + theme(legend.position='hidden'),
grid.arrange(p2a$plot + theme(legend.position='hidden'),p2a$table, layout_matrix = rbind(c(1), c(1), c(2))),
p2b + theme(legend.position='hidden'),
grid.arrange(p3a$plot + theme(legend.position='hidden'),p3a$table, layout_matrix = rbind(c(1), c(1), c(2))),
p3b + theme(legend.position='hidden'),
ncol = 2)


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
	ggsurvplot(fit1, data = dataset, pval = T, risk.table= T,pval.size =6, title= datasettitle, ggtheme=custom_theme())
}

make_tmb_group_dataset <- function(dataset) {
	dataset <- dataset %>% filter(!is.na(tmb)) 
	tmb_med <- dataset %>% select(tmb) %>% as.matrix() %>% as.vector() %>% median()
	dataset <- dataset %>% mutate(Known_TMB = ifelse(tmb >= tmb_med, "High", "Low"))
	return(dataset)
}

my_tmb_survival_plot <- function(dataset, datasettitle) {
	fit1 <- survfit(Surv(time = time , event = status )~ Known_TMB , data = dataset)
	ggsurvplot(fit1, data = dataset, pval = T, pval.size =6, risk.table= T, title= datasettitle, ggtheme=custom_theme())
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
	data$response <- factor(data$response, levels = c("R", "NR"))

	ht <- sum(data[data$Novel_Signature == "High",]$n)
	lt <- sum(data[data$Novel_Signature == "Low",]$n)
	data$tt <- c(rep(ht, nrow(data[data$Novel_Signature == "High",])), rep(lt,nrow(data[data$Novel_Signature == "Low",])))
	data <- data %>% mutate(percent = round(n/tt*100,1)) 

	dt <- dcast(data, Novel_Signature ~ response.binary, value.var="n")
	change_na2zero <- function(x) {ifelse(is.na(x), 0, as.double(x))}
	apply(dt[,2:5], 2, change_na2zero) %>% fisher.test() -> dtt

	# make plot
	p <- ggplot(data, aes(fill=response.binary, y=n, x=Novel_Signature)) +
    geom_bar(position="fill", stat="identity") +
		scale_x_discrete(labels=c(paste0("High\n(n=", ht,")"), paste0("Low\n(n=", lt,")")))

	p <- p + scale_fill_manual(values = c("orange", "lightskyblue3") +
	theme_classic() + xlab("") + ylab("Percentage of patients")
	p <- p + geom_text(aes(label=paste0(n,"(",percent,"%)")), position = position_fill(vjust = 0.5))+
	guides(fill=guide_legend(title="Responder"))+theme(axis.text.x = element_text(size = 12))

	p <- p + scale_y_continuous(breaks = c(0, 0.5, 1,1.5))+
                             geom_signif(annotations = paste0("pvalue=", round(dtt$p.value,2)),
                                         y_position = 1.03)
	geom_text(x =1.8, y =1.03, label = paste0("pvalue=", round(dtt$p.value,2)), hjust = 1, size = 4)
	
	return(p)
}


### 02. Load data
load(paste0(indir,"clinical_gse176307.RData"))
load(paste0(indir,"clinical_imvigor210core.RData"))
load(paste0(indir,"clinical_ucgenome.RData"))



### 03. Survival relevance and ICB response of novel bladder signature in GSE176307
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



### 04. Survival relevance and ICB response of novel bladder signature in IMvigor210


### 05. Survival relevance and ICB response of novel bladder signature in UC-GENOME


### 06. Survival relevance and ICB response of known TMB in GSE176307


### 07. Survival relevance and ICB response of known TMB in IMvigor210


### 08. Survival relevance and ICB response of known TMB in UC-GENOME

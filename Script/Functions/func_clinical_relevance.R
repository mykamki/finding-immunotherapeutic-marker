
tpm <- function(counts, lengths) {
	rpk <- counts/(lengths/1000) # 1. rpk
	coef <- sum(rpk)/1e6 # 2. coef(per million scaling factor)
	rpk/coef # 3. tpm
}

zscore_transform <- function(x) {
  (x-mean(x))/sd(x)
}

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}

my_geneset_survival_plot <- function(dataset, datasettitle) {
	fit1 <- survfit(Surv(time = time , event = status )~ Novel_Signature , data = dataset)
	ggsurvplot(fit1, data = dataset, pval = T, pval.size =6, risk.table= T, 
                   title= datasettitle, legend = "right",
                   ggtheme=custom_theme(), tables.height = 0.1,
                   tables.theme = theme_cleantable(), fontsize = 3, 
                   risk.table.y.text.col = T, # colour risk table text annotations.
                   risk.table.y.text = F, # show bars instead of names in text annotations
                   #legend.labs = expression("Signature[High]", "Signature[Low]")  # Change legend labels
                )
}

make_tmb_group_dataset <- function(dataset) {
	dataset <- dataset %>% filter(!is.na(tmb)) 
	tmb_med <- dataset %>% select(tmb) %>% as.matrix() %>% as.vector() %>% median()
	dataset <- dataset %>% mutate(Known_TMB = ifelse(tmb >= tmb_med, "High", 
							 ifelse(tmb < tmb_med, "Low", NA)))
	return(dataset)
}

my_tmb_survival_plot <- function(dataset, datasettitle) {
	fit1 <- survfit(Surv(time = time , event = status )~ Known_TMB , data = dataset)
	ggsurvplot(fit1, data = dataset, pval = T, pval.size =6, risk.table= T, 
                   title= datasettitle, legend = "right",
                   ggtheme=custom_theme(), tables.height = 0.1,
                   tables.theme = theme_cleantable(), fontsize = 3, 
                   risk.table.y.text.col = T, # colour risk table text annotations.
                   risk.table.y.text = F, # show bars instead of names in text annotations
                   #legend.labs = expression("Signature[High]", "Signature[Low]")  # Change legend labels
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
                     axis.text.y = element_blank(),legend.position="top",
                     plot.margin = margin(0, 1, 0, 1, "cm")) +
                geom_text(aes(label=paste0(n," (",percent,"%)")), position = position_stack(vjust = 0.5)) +
                guides(fill=guide_legend(title="Responder")) + 
                theme(axis.text.x = element_text(size = 12))
	
	return(p)
}

bar_plot_by_tmb <- function(dataset) {
	data <- dataset %>% 
                filter(!is.na(response.binary)) %>% 
		filter(!is.na(Known_TMB)) %>% 
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
                     axis.text.y = element_blank(), legend.position="top",
                     plot.margin = margin(0, 1, 0, 1, "cm")) +
                geom_text(aes(label=paste0(n," (",percent,"%)")), position = position_stack(vjust = 0.5)) +
                guides(fill=guide_legend(title="Responder")) + 
                theme(axis.text.x = element_text(size = 12))
	
	return(p)
}

extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

			

my_anno_plot <- function(mylabel, fontsize = 6) {
data <- as.data.frame(matrix(seq(1,20), nrow =4))
# Create the plot with modified axis settings
ggplot(data, aes(x = V1, y = V2)) + xlab("") +
  theme(
    axis.line.x.top = element_line(colour = "white"), 
    axis.line.y.left = element_line(colour = "white"), 
    axis.line.y.right = element_line(colour = "white"), 
    axis.line.x.bottom = element_line(colour = "black"), 
    axis.text.y = element_blank(),  # Hide y-axis labelsgeomgegg
    axis.ticks.y = element_blank(),  # Hide y-axis ticks
    axis.title.y = element_blank(),  # Hide y-axis title
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white")
  ) +
  scale_x_continuous(breaks = NULL) + # Hide x-axis ticks
geom_text(aes(x = 2, y = 4.5), label = mylabel, size=fontsize)
}


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






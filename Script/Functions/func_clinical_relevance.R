
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
                   title= datasettitle, 
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
                   title= datasettitle, 
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

			

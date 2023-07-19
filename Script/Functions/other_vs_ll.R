other_vs_ll <- function() 
  {
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

### Survival relevance and ICB response of combination between novel bladder signature and known TMB in GSE176307
p1a <- my_combination_survival_plot(clinical_gse176307, "GSE176307")
p1b <- bar_plot_by_combination(clinical_gse176307)


### Survival relevance and ICB response of ncombination between novel bladder signature and known TMB in IMvigor210
p2a <- my_combination_survival_plot(clinical_imvigor210core, "IMvigor210")
p2b <- bar_plot_by_combination(clinical_imvigor210core)


### Survival relevance and ICB response of combination between novel bladder signature and known TMB in UC-GENOME
p3a <- my_combination_survival_plot(clinical_ucgenome, "UC-GENOME")
p3b <- bar_plot_by_combination(clinical_ucgenome)


###  Plotting
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
			

p <- grid.arrange(pA, shared_legend, heights = c(8, 1))

return(p)
  }

hh_vs_other <- function(clinical_gse176307, clinical_imvigor210core, clinical_ucgenome) 
{
  ### Survival relevance and ICB response of combination between novel bladder signature and known TMB in GSE176307
  p1a <- my_combination_survival_plot(clinical_gse176307, "GSE176307")
  p1b <- bar_plot_by_combination(clinical_gse176307)


  ### Survival relevance and ICB response of ncombination between novel bladder signature and known TMB in IMvigor210
  p2a <- my_combination_survival_plot(clinical_imvigor210core, "IMvigor210")
  p2b <- bar_plot_by_combination(clinical_imvigor210core)


  ### Survival relevance and ICB response of combination between novel bladder signature and known TMB in UC-GENOME
  p3a <- my_combination_survival_plot(clinical_ucgenome, "UC-GENOME")
  p3b <- bar_plot_by_combination(clinical_ucgenome)

  ### Plotting
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

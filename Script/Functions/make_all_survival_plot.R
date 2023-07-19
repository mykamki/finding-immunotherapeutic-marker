make_all_survival_plot <- function(clinical_gse176307, clinical_imvigor210core, clinical_ucgenome) 
{
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

  comper <- grid.arrange(qa, qb, qc, shared_legend, heights = c(2,2,2,1))
  return(comper)
  }

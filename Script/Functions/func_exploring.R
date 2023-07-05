make_cluster_n <- function(which_scissor_cell) {
	df <- as.data.frame(data1$seurat_clusters[names(data1$seurat_clusters) %in% which_scissor_cell])
	colnames(df) <- "clu"
	dfs <- df %>% group_by(clu) %>% summarize(n = n()) %>% arrange(-n)
	dfs$name <- paste0("cluster", dfs$clu)
	dfs <- dfs[,-1]
	return(dfs)
	}

make_celltype_n <- function(which_scissor_cell) {
	data3 <- subset(data1, subset = CellName %in% which_scissor_cell)
	df2 <- as.data.frame(table(Idents(data3)))
	colnames(df2) <- c("CellType","No")
	df2 <- df2 %>% arrange(-No)
	df2$name <- gsub("_", " ", df2$CellType, fixed = T)
	return(df2)
	}


test_de_cutoff <- function(de, threshold) {
catego_updown <- function(de.res) {
  de.res <- de.res %>%
    mutate(
    expression = case_when(avg_log2FC >= threshold & p_val_adj <= 0.05 ~ "Up-regulated",
                           avg_log2FC <= -threshold & p_val_adj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
    )
  return(de.res)
}
de <- catego_updown(de)
}

get_top10 <- function(de.res) {
  top <- 20
  top_genes <- bind_rows(
  de.res %>% 
    filter(expression == 'Up-regulated') %>% 
    arrange(desc(abs(avg_log2FC)), p_val_adj) %>% 
    head(top),
  de.res %>% 
    filter(expression == 'Down-regulated') %>% 
    arrange(desc(abs(avg_log2FC)), p_val_adj) %>% 
    head(top)
  )    
  top_genes <- top_genes %>% mutate (Genes = rownames(top_genes))
  return(top_genes)
}

my_volcano <- function(de.res) {
  top_genes <- get_top10(de.res)
  sigupn <- nrow(de.res[de.res$expression == "Up-regulated",])
  sigdownn <- nrow(de.res[de.res$expression == "Down-regulated",])

  p <- ggplot(de.res, aes(avg_log2FC, -log(p_val_adj,10))) +
    geom_point(aes(color = expression), size = 3, alpha = 0.7) +
		scale_x_continuous(limits = c(-10, 10))+
	  scale_y_continuous(limits = c(0, 200))+
    xlab(expression("log"[2]*"FC")) +
    ylab(expression("-log"[10]*"FDR")) +
    scale_color_manual(values = c( "dodgerblue3", "gray50", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    geom_hline(yintercept = -log(0.05,10)) +
    geom_vline(xintercept = 1.5, linetype = "dotted") +
    geom_vline(xintercept = -1.5, linetype = "dotted") +
    theme_bw()+
      theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_blank(),
				panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

  p <- p + geom_label_repel(data = top_genes,
                   mapping = aes(avg_log2FC, -log(p_val_adj,10), label = Genes),
                   size = 2, max.overlaps= 30)

  #xpos <- range(de.res$avg_log2FC)*1/2
  #ypos <- range(-log(de.res$p_val_adj,10))*0.92

	#xpos <- c(-10,10)*1/2
	#ypos <- c(-1,200)*5/9

  #p <- p + geom_text(x = xpos[1], y = ypos[2],
	#           label = paste0("Down (", sigdownn, ")"), col = "dodgerblue3", size = 4) +
  #         geom_text(x = xpos[2], y = ypos[2],
  #           label = paste0("Up (", sigupn, ")"), col = "firebrick3", size = 4)

	p <- p + theme(legend.position= c(0.2,0.85))
  return(p)
}



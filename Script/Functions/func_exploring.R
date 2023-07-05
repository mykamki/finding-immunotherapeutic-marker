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

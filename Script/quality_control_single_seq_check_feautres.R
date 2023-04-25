library(Seurat)
indir <- ""
pid <- ""



### 01. Load data 
sc_files <- list.files(path = location, pattern = paste0(pid, "_gene_cell_exprs_table"))
scdata <- fread(paste0(location, sc_files[1]))


#### 2. remove duplicate ####

sc_dataset <- process_scdataset(scdata)


process_scdataset <- function(dataset) {
	a <- dataset%>% group_by(Symbol) %>% summarize(n= n()) %>% filter(n>1)
	idx <- list()
	for (n in 1:length(a$Symbol)) {
		b <- which(a$Symbol[n] == dataset$Symbol)[2]
		idx <- append(b, idx)
	}
	idx <- as.double(flatten_chr(idx))
	b <- dataset[!idx,]
  c <- as.matrix(b[,-c(1,2)])
	rownames(c) <- b$Symbol
   colnames(c) <- colnames(b)[-c(1,2)]

	return(c)
}

number <- rev(strsplit(rev(strsplit(args[2], "/")[[1]])[1], "sc")[[1]])[1]

cat("#####################################################\n")
cat("## The number of gene of sc_dataset : ", dim(sc_dataset)[1], " ##\n")
cat("## The number of cell of sc_dataset : ", dim(sc_dataset)[2], " ##\n")
cat("#####################################################\n")
save(sc_dataset, file = paste0(args[2], "sc_dataset_count.RData"))


#### 3. Prepare for seurat format ####
sc_dataset <- Seurat_preprocessing(sc_dataset, verbose = F)


#### 4. Save data ####
save(sc_dataset,  file = paste0(args[2], "sc_dataset.RData"))

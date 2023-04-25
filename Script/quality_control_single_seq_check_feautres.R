library(Seurat)
indir <- "../Data/"
pid <- ""
outdir <- "../Output/"


### 01. Load data 
sc_files <- list.files(path = location, pattern = paste0(pid, "_gene_cell_exprs_table"))
scdata <- fread(paste0(location, sc_files[1]))



### 02. processing data
sc_dataset <- process_scdataset(scdata)
cat("## The number of gene of sc_dataset : ", dim(sc_dataset)[1], " ##\n")
cat("## The number of cell of sc_dataset : ", dim(sc_dataset)[2], " ##\n")



### 03. make Seurat object
data <- CreateSeuratObject(counts = sc_dataset, project = args[2], min.cells = 3,  min.features= 200)


### 04. make mitochondrial gene
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-") 


### 05. plot the features
# Visualize QC metrics as a violin plot
pdf(file = paste0(outdir, "qcplot1.pdf"))
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
  
#FeatureScatter is typically used to visualize feature-feature relationships, but can be used
#for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(file = paste0(outdir, "qcplot2.pdf"))
plot1 + plot2 
dev.off()


### 06. save the results
save(sc_dataset, file = paste0(outdir, "sc_dataset_count.RData"))




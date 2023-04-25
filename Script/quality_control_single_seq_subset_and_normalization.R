
library(Seurat)


#### 01. Load data ####
load(file = paste0(args[1], "sc_dataset_count.RData"))
data <- CreateSeuratObject(counts = sc_dataset, project = args[2], min.cells = 3,  min.features= 200) # make Seurat object
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-") # make mitochondrial gene


#### 02. Subset data ####
min.nfeature <- as.double(args[3])
max.nfeature <- as.double(args[4])
max.mt <- as.double(args[5])

data1 <- subset(data, subset = nFeature_RNA > min.nfeature & nFeature_RNA < max.nfeature & percent.mt < max.mt)  # min.cells = 400, min.features = 0


cat("################################################", "\n")
cat("## The number of cell before QC : ", length(data$orig.ident), " ##", "\n")
cat("## The number of cell after QC : ", length(data1$orig.ident), " ##", "\n")
cat("################################################", "\n")



#### 03. Normalizing the data ####
data1 <- NormalizeData(data1, normalization.method = "LogNormalize", scale.factor = 10000)
data1 <- FindVariableFeatures(data1, selection.method = "vst", nfeatures = 2000)



#### 04. Scaling the data ####
all.genes <- rownames(data1)
data1 <- ScaleData(data1, features = all.genes)



#### 05. Perform linear dimensional reduction ####
data1 <- RunPCA(data1, features = VariableFeatures(object = data1))
pdf(file = paste0(args[6], "elbowplot.pdf"))
ElbowPlot(data1, ndims= 50)
dev.off()

save(data1, file = paste0(args[6], "data1.RData"))

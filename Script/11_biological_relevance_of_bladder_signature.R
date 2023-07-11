
### 01. Setting
library(ComplexHeatmap)
library(circlize)
library(dplyr)


### 02. Load pathway gene list
path.genes <- read.csv("/storage/home/mina_20/ANALYSIS/singlecell/analysis_script/15_IMvgor210/path.genes.csv", header = T)
path.genes$group <- gsub("Immune_checkpoint", "CheckPoint", path.genes$group)
path.genes$group <- gsub("Cell_cycle", "cellcycle", path.genes$group)


### 03. make annotation data
cl <- clinical_gse176307[,c("ID", "Novel_Signature", "Known_TMB", "response.binary")]
cl <- cl[order(cl$Novel_Signature, cl$Known_TMB),]

ann_colors = list(Novel_Signature = c("High" = "coral2", "Low" = "dodgerblue1"),
								Known_TMB =  c("High" = "coral2", "Low" = "dodgerblue1"),
               response.binary = c("NR" = "lightskyblue3", "R" = "orange",  "NA" = "gray60")
							)
anno <- cl[,-1]

ha1 = HeatmapAnnotation(df = anno,
    col = ann_colors
		)


### 04. make heatmap matrix
myz <- function(x) {
	(x-mean(x))/sd(x)
}

mat <- log_dataset_gse176307[rownames(log_dataset_gse176307) %in% path.genes$gene,]
mat <- apply(mat,1, myz)
mat <- t(mat)

if (identical(cl$ID, colnames(mat))) {
  print("ok")
  } else {
  mat <- mat[,cl$ID]
}
identical(cl$ID, colnames(mat)) # check again


### 05. Draw heatmap
for (n in 1:length(unique(path.genes$group))) {
    path <- unique(path.genes$group)[n]
    submat <- mat[rownames(mat) %in% path.genes[path.genes$group == path,]$gene,]
    ht <- Heatmap(submat, cluster_columns = FALSE, column_labels= FALSE, show_column_names = FALSE, row_title = path, show_heatmap_legend = FALSE)
    assign(paste0("ht",n), ht)
  }

path <- unique(path.genes$group)[1]
submat <- mat[rownames(mat) %in% path.genes[path.genes$group == path,]$gene,]
ht1 <- Heatmap(submat, cluster_columns = FALSE, column_labels= FALSE, show_column_names = FALSE, row_title = path, name = "Gene expression")


ht_list = ha1 %v% ht1 %v% ht3 %v% ht2 %v% ht4 %v% ht5 %v% ht6 %v% ht9 %v% ht7 %v% ht8 %v% ht10 %v% ht11 
png(paste0(outdir, "heatmap.png"),width = 500, height = 850); draw(ht_list, ht_gap = unit(0.3, "cm")) ;dev.off()


### 06. T-test for pathway
for (n in 1:length(unique(path.genes$group))) {
    path <- unique(path.genes$group)[n]
    submat <- mat[rownames(mat) %in% path.genes[path.genes$group == path,]$gene,]

    high_novel_gene_signature <- apply(submat[,colnames(submat) %in% cl[cl$Novel_Signature == "High",]$ID], 2, mean)
    low_novel_gene_signature <- apply(submat[,colnames(submat) %in% cl[cl$Novel_Signature == "Low",]$ID], 2, mean)
   
    high_known_tmb <- apply(submat[,colnames(submat) %in% cl[cl$Known_TMB == "High",]$ID], 2, mean)
    low_known_tmb <- apply(submat[,colnames(submat) %in% cl[cl$Known_TMB == "Low",]$ID], 2, mean)

    #print(path)
    #t.test(high_novel_gene_signature, low_novel_gene_signature) %>% print()
    #t.test(high_known_tmb, low_known_tmb) %>% print()

    df <- data.frame(value = c(high_novel_gene_signature, low_novel_gene_signature, high_known_tmb, low_known_tmb),
                 marker = c(rep("Novel gene signature", length(high_novel_gene_signature)+length(low_novel_gene_signature)),
                            rep("Known TMB", length(high_known_tmb)+length(low_known_tmb))),
                  group = c(rep("High", length(high_novel_gene_signature)), 
                            rep("Low", length(low_novel_gene_signature)),
                            rep("High", length(high_known_tmb)), 
                            rep("Low", length(low_known_tmb)))
                 )
    df$marker <- factor(df$marker, levels = c("Novel gene signature", "Known TMB")) 
    assign(paste0("bx",n), genepath_boxplot(df))  
  
  }

genepath_boxplot <- function(data) {
  ggplot(data, aes(x=group, y=value, fill=group)) +
  geom_boxplot() +
  facet_wrap(~marker) +
  scale_colour_manual(values = c("coral2", "dodgerblue1")) +
  theme_bw() + 
	scale_fill_manual(values = c("coral2", "dodgerblue1")) +
	xlab("") + ylab(paste0("The average of","\n", "genes expression in pathway")) + ggtitle(path) +
	geom_signif(comparisons = list(c("High", "Low")), map_signif_level = TRUE, step_increase = 0.1) +
	theme_half_open(12) +
  theme(plot.title = element_text(hjust = 0.5),
		legend.position="rigth",
    legend.title = element_blank(),
		strip.background = element_rect(fill="white",colour="white"))
  }

png(paste0(outdir, "test.png")); genepath_boxplot(df) ;dev.off()

library(cowplot)
library(grid)
hts <- draw(ht_list, ht_gap = unit(0.3, "cm"))
htsg <- grid.grabExpr(ht_list, ht_gap = unit(0.3, "cm"))
bxs <- grid.arrange(bx2 +theme(legend.position="none"), 
                    bx3 +theme(legend.position="none"), 
                    bx4 +theme(legend.position="none"),
                    bx11 +theme(legend.position="none"), ncol =2)

png(paste0(outdir, "bx.png"), width = 300, height = 1200)
grid.arrange(bx2 +theme(legend.position="none"), 
                    bx3 +theme(legend.position="none"), 
                    bx4 +theme(legend.position="none"),
                    bx11 +theme(legend.position="none"), ncol =1)
dev.off()

png(paste0(outdir, "supfig7.png"), width = 800, height = 900)
grid.arrange(bx1 +theme(legend.position="none"), 
            bx5 +theme(legend.position="none"), 
            bx6 +theme(legend.position="none"),
              bx7 +theme(legend.position="none"),
              bx8 +theme(legend.position="none"),
              bx9 +theme(legend.position="none"),
              bx10 +theme(legend.position="none"), ncol =3)
dev.off()

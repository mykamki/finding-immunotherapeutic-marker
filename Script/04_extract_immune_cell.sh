# annotation marker to integrated dataset
# outdir is "cell_type_ann_result"
Rscript make_violin_plot_using_cell_markers.R ../Output/INTEGRATION/

# mannually
vim ../Output/INTEGRATION//cell_type_ann_result/cluster_ann.R 

# re-assign cell type using cell-marker genes
mkdir ~/ANALYSIS/singlecell/Rdata/integ_3000_mt/dim30/anned_data1 
Rscript 3_MARKER/marker02_reassign_celltype.R ~/ANALYSIS/singlecell/Rdata/integ_3000_mt/dim30/ ~/ANALYSIS/singlecell/Rdata/integ_3000_mt/dim30/cell_type_ann_result/cluster_ann.R ~/ANALYSIS/singlecell/Rdata/integ_3000_mt/dim30/anned_data1/ 

# split data for purpose
mkdir ~/ANALYSIS/singlecell/Rdata/integ_3000_mt/dim30/rm_epi
Rscript 4_SUBSET/subset01_rm_epithelial.R ~/ANALYSIS/singlecell/Rdata/integ_3000_mt/dim30/anned_data1/ ~/ANALYSIS/singlecell/Rdata/integ_3000_mt/dim30/rm_epi/  

####################################################
#### Before removing  Epithelial_cell  :  36300 ####
##### After removing  Epithelial_cell  :  1517 #####
####################################################

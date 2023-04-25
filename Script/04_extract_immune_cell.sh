# annotation marker to integrated dataset
# outdir is "cell_type_ann_result"
Rscript make_violin_plot_using_cell_markers.R ../Output/INTEGRATION/

# mannually
vim ../Output/INTEGRATION/cell_type_ann_result/cluster_ann.R 

# re-assign cell type using cell-marker genes and extract just immune cells
mkdir ../Output/INTEGRATION/IMMUNE/
Rscript remove_specific_cell_type.R ../Output/INTEGRATION/../Output/INTEGRATION/cell_type_ann_result/cluster_ann.R Epithelial_cell ../Output/INTEGRATION/IMMUNE/

####################################################
#### Before removing  Epithelial_cell  :  36300 ####
##### After removing  Epithelial_cell  :  1517 #####
####################################################


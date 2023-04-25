# Integration using SCT
# Integration goals
# 1. Create an ‘integrated’ data assay for downstream analysis
# 2. Identify cell types that are present in both datasets
# 3. Obtain cell type markers that are conserved in both control and stimulated cells
# 4. Compare the datasets to find cell-type specific responses to stimulation


# 01. Integration 
Rscript integration_single_seq_by_SCT.R ../Output/ ../Output/INTEGRATION/


# 02. Check optimal resolution for clustering
Rscript quality_control_single_seq_find_optimal_resolution.R ../Output/INTEGRATION/ 30


# 03. Do clustering
Rscript 2_INTEGRATION/integ05_clustering.R ~/ANALYSIS/singlecell/Rdata/integ_3000_mt/ 30 1 ~/ANALYSIS/singlecell/Rdata/integ_3000_mt/dim30/



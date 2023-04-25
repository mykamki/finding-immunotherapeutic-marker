# preprocessing sc data

# 1. make dataset and check nFeature for determine minimal cell counts
for n in {1..7};do Rscript quality_control_single_seq_check_feautres.R ../Data/ SC${n} ../Output/SC${n]/; done


# 2. input hyper-parameter for determine minimal cell counts
Rscript quality_control_single_seq_subset_and_normalization.R ../Output/SC1/ BC1 200 7500 10 
Rscript quality_control_single_seq_subset_and_normalization.R ../Output/SC2/ BC2 1000 6000 10 
Rscript quality_control_single_seq_subset_and_normalization.R ../Output/SC3/ BC3 500 6000 10 
Rscript quality_control_single_seq_subset_and_normalization.R ../Output/SC4/ BC4 700 3000 10 
Rscript quality_control_single_seq_subset_and_normalization.R ../Output/SC5/ BC5 500 4000 10 
Rscript quality_control_single_seq_subset_and_normalization.R ../Output/SC6/ BC6 500 5000 7.5
Rscript quality_control_single_seq_subset_and_normalization.R ../Output/SC7/ BC7 400 4000 10 


# 3. input hyper-parameter for pca and find optimal resolution
for n in {1..7};do Rscript quality_control_single_seq_find_optimal_resolution.R ../Output/SC${n]/ 10 ; done
Rscript quality_control_single_seq_find_optimal_resolution.R ~/ANALYSIS/singlecell/Rdata/sc5/ 15 


# 4. input hyper-parameter for clustering
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc1/ 10 1.6 ~/ANALYSIS/singlecell/Rdata/sc1/  
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc2/ 10 0.5 ~/ANALYSIS/singlecell/Rdata/sc2/ 
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc3/ 10 1.6 ~/ANALYSIS/singlecell/Rdata/sc3/ 
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc4/ 10 0.5 ~/ANALYSIS/singlecell/Rdata/sc4/
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc5/ 15 2.0 ~/ANALYSIS/singlecell/Rdata/sc5/ 
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc6/ 10 4.0 ~/ANALYSIS/singlecell/Rdata/sc6/ 
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc7/ 10 2.0 ~/ANALYSIS/singlecell/Rdata/sc7/

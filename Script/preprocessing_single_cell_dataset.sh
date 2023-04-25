# preprocessing sc data

# 1. make count dataset
for n in {1..7};do Rscript quality_control_single_seq_check_feautres.R ../Data/ SC${n} ../Output/SC${n]/; done

# 1. make dataset and check nFeature for determine minimal cell counts
for n in {1..7};do Rscript 1_QC/qc01_check_nFeature.R ~/ANALYSIS/singlecell/Rdata/sc$n/ BC$n ~/ANALYSIS/singlecell/Rdata/sc$n/; done

# 2. input hyper-parameter for determine minimal cell counts
Rscript 1_QC/qc02_subset_and_normalization.R ~/ANALYSIS/singlecell/Rdata/sc1/ BC1 200 7500 10 ~/ANALYSIS/singlecell/Rdata/sc1/
Rscript 1_QC/qc02_subset_and_normalization.R ~/ANALYSIS/singlecell/Rdata/sc2/ BC2 1000 6000 10 ~/ANALYSIS/singlecell/Rdata/sc2/
Rscript 1_QC/qc02_subset_and_normalization.R ~/ANALYSIS/singlecell/Rdata/sc3/ BC3 500 6000 10 ~/ANALYSIS/singlecell/Rdata/sc3/
Rscript 1_QC/qc02_subset_and_normalization.R ~/ANALYSIS/singlecell/Rdata/sc4/ BC4 700 3000 10 ~/ANALYSIS/singlecell/Rdata/sc4/
Rscript 1_QC/qc02_subset_and_normalization.R ~/ANALYSIS/singlecell/Rdata/sc5/ BC5 500 4000 10 ~/ANALYSIS/singlecell/Rdata/sc5/
Rscript 1_QC/qc02_subset_and_normalization.R ~/ANALYSIS/singlecell/Rdata/sc6/ BC6 500 5000 7.5 ~/ANALYSIS/singlecell/Rdata/sc6/
Rscript 1_QC/qc02_subset_and_normalization.R ~/ANALYSIS/singlecell/Rdata/sc7/ BC7 400 4000 10 ~/ANALYSIS/singlecell/Rdata/sc7/

# 3. input hyper-parameter for pca and find optimal resolution
for n in {1..7};do Rscript 1_QC/qc03_find_optimal_resolution.R ~/ANALYSIS/singlecell/Rdata/sc$n/ 10 ~/ANALYSIS/singlecell/Rdata/sc$n/; done
Rscript 1_QC/qc03_find_optimal_resolution.R ~/ANALYSIS/singlecell/Rdata/sc5/ 15 ~/ANALYSIS/singlecell/Rdata/sc5/

# 4. input hyper-parameter for clustering
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc1/ 10 1.6 ~/ANALYSIS/singlecell/Rdata/sc1/  
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc2/ 10 0.5 ~/ANALYSIS/singlecell/Rdata/sc2/ 
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc3/ 10 1.6 ~/ANALYSIS/singlecell/Rdata/sc3/ 
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc4/ 10 0.5 ~/ANALYSIS/singlecell/Rdata/sc4/
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc5/ 15 2.0 ~/ANALYSIS/singlecell/Rdata/sc5/ 
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc6/ 10 4.0 ~/ANALYSIS/singlecell/Rdata/sc6/ 
Rscript 1_QC/qc04_clustering.R ~/ANALYSIS/singlecell/Rdata/sc7/ 10 2.0 ~/ANALYSIS/singlecell/Rdata/sc7/

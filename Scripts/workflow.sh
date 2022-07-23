#!/bin/bash

DIR="/data/public/adesous1/GeneCorrelation" # work directory
EXTDIR="/data/public/adesous1/GeneCorrelation/External" # where to fetch external data & scripts (network raw files, GTEx, Tabula Sapiens, etc)
cd $DIR


## 0. Preprocessing -----------------------------------------------------------

for NET in "stabsel" "stabsel_pcclasso" "stabsel_randomized"
do
  echo $NET
  # Grab "raw" networks from external location, read and save appropriately
  ./Scripts/fetch_networks.sh ${EXTDIR} ${NET} "0_Preprocessing"
  # Transform networks (remove residual edges, reduce to largest connected
  # component)
  ./Scripts/run_preprocessing.sh ${NET} "0_Preprocessing"
done


## 1. Topology analysis -------------------------------------------------------

for NET in "stabsel" "stabsel_filtered" "stabsel_filtered_largestCC"  "stabsel_filtered_trans_largestCC" "stabsel_pcclasso" "stabsel_pcclasso_filtered" "stabsel_pcclasso_filtered_largestCC"
do
  echo $NET
  ./Scripts/run_topology_analysis.sh ${NET} "1_Topology"
done


## 2. Clustering --------------------------------------------------------------

for NET in "stabsel_filtered_largestCC"
do
  echo $NET
  ./Scripts/run_clustering_analysis.sh ${NET} "2_Clustering"
done


## 3. Expression data preparation ---------------------------------------------

for DATA in GTEx TS
do
  echo $DATA
  ./Scripts/run_${DATA}dataprep.sh "3_${DATA}DataPrep"
done


## 4. Module activity ---------------------------------------------------------

for NET in "stabsel_filtered_largestCC"
do
  echo $NET
  ./Scripts/run_activity_analysis.sh ${NET} "4_ModuleActivity"
done


# ## Module activity ------------------------------------------------------------
# 
# 
# # Weight LFCs based on module's local topology
# Scripts/${WDIR}/Compute_LocallyWeightedActivities.R
# 
# # Smooth GTEx LFCs through all modules
# Scripts/${WDIR}/SmoothLFCs_MM.R
# 
# # Average smoothed values across gene modules
# Scripts/${WDIR}/Average_SmoothedLFCs.R
# 
# # Plot obtained relative activities
# Scripts/${WDIR}/Plot_RelativeActivities.R
# 
# 
# # ## Predictability analysis
# # 
# # # First, predict expression values using the networks in GTEx subsets:
# # predict_expression.R
# # 
# # # Then, compute correlation between predicted and observed values:
# # compute_network_correlations.R
# # 
# # # This is done per tissue and for the cross-tissue (all subsets of the same size).
# # 
# # # Then, analyse overlap of poorly predicted genes across tissues:
# # compute_plot_poorly_predicted.R
# # 
# # # Characterize the poorly predicted genes based on dataset- and network-specific features:
# # characterize_poorly_predicted.R
# # 
# # # Compare predictability with and without partial correlations:
# # plot_predictability_comparison_GTEx_subsets.R
# # plot_predictability_stabsel_on_pcclasso.R
# 


## Predictability changes with age

./Scripts/run_predictability_analysis.sh

# # Compute difference between predicted by network and observed, for each gene in each sample:
# compute_network_delta.R
# 
# # Limma for both expression and predictability changes with age:
# compute_age_limma_GTEx.R
# 
# # Plot expression vs predictability:
# plot_age_cordeltas_logFCs_tissues.R
# 
# # GO enrichment in genes with significant predictabilty changes with age:
# expression_predictability_FCs.R

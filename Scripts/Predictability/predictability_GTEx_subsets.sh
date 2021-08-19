#!/bin/bash

# Predict expression in GTEx subsets divided by tissue, age and both
# Using both full network and network with no indirect effects


for NET in stabsel stabsel_pcclasso
do
  for TYPE in Age Tissue AgeTissue
  do
    # predict expression on GTEx subsets
    Rscript --vanilla predict_expression.R $NET $TYPE "T"
    Rscript --vanilla predict_expression.R $NET $TYPE "F"
    
    # compute net correlations based on previous predictions
    Rscript --vanilla compute_network_correlations.R $NET $TYPE "T"
    Rscript --vanilla compute_network_correlations.R $NET $TYPE "F"
  done
done


# plot comparisons
for TYPE in Age Tissue AgeTissue
do
  Rscript --vanilla plot_predictability_comparison_GTEx_subsets.R $TYPE "T"
done
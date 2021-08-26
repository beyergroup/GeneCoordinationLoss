#!/bin/bash

for NET in stabsel stabsel_pcclasso
do
  # compute difference between network predictions and original data
  Rscript --vanilla compute_network_delta.R $NET
  
  # limma for age-related changes in predictability (DP)
  Rscript --vanilla compute_age_limma_GTEx.R $NET "delta"
  
  # limma for age-related DE
  Rscript --vanilla compute_age_limma_GTEx.R $NET "expression"
  
  # plot DE vs DP for all tissues, excluding poorly predicted genes
  Rscript --vanilla plot_age_cordeltas_logFCs_tissues.R $NET "T"
  
  # GO enrichmemt for significantly DP genes
  Rscript --vanilla GO_enrichment_DP.R $NET "T"
done

#!/bin/bash

NET=$1
WDIR=$2
echo "Mode = ${NET}"

mkdir Outputs/${WDIR}
mkdir Plots/${WDIR}
mkdir Logs/${WDIR}

SUBDIR="PoorlyPredicted"
mkdir Outputs/${WDIR}/${SUBDIR}
mkdir Plots/${WDIR}/${SUBDIR}
mkdir Logs/${WDIR}/${SUBDIR}

NET="stabsel_filtered_trans_largestCC"

# Compute correlation thresholds
Rscript --vanilla Scripts/${WDIR}/${SUBDIR}/Compute_WellPredicted.R "sampled_centered_correlations" Outputs/${WDIR}/${SUBDIR} Outputs/${WDIR}/${SUBDIR}/Randomized "Spearman" &> Logs/${WDIR}/${SUBDIR}/Compute_WellPredicted_sampled_centered_Spearman_correlations_$(date +%F)--$(date +%H)":"$(date +%M).out

# Characterize poorly predicted genes
Rscript --vanilla Scripts/${WDIR}/${SUBDIR}/Characterize_PoorlyPredicted.R <!GTEX DIR!> Outputs/${WDIR}/${SUBDIR} "sampled_centered_correlations_spearman" &> Logs/${WDIR}/${SUBDIR}/Characterize_PoorlyPredicted_sampled_centered_Spearman_correlations_$(date +%F)--$(date +%H)":"$(date +%M).out

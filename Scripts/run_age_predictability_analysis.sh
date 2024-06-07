#!/bin/bash

NET=$1
WDIR=$2
echo "Mode = ${NET}"

mkdir Outputs/${WDIR}
mkdir Plots/${WDIR}
mkdir Logs/${WDIR}

SUBDIR="Age"
mkdir Outputs/${WDIR}/${SUBDIR}
mkdir Plots/${WDIR}/${SUBDIR}
mkdir Logs/${WDIR}/${SUBDIR}

NET="stabsel_filtered_trans_largestCC"

# True Network predictions per age group and tissue
Rscript --vanilla Scripts/${WDIR}/Predict_Expression.R ${NET} "sampled_centered" "Outputs/3_GTExDataPrep/Subset_Data" Outputs/${WDIR}/${SUBDIR} &> Logs/${WDIR}/${SUBDIR}/Predict_Expression_${NET}_Age_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Predict_Expression.R ${NET} "maxsampled_centered" "Outputs/3_GTExDataPrep/Subset_Data/Max_Subset" Outputs/${WDIR}/${SUBDIR}/Age_MaxSubset &> Logs/${WDIR}/${SUBDIR}/Predict_Expression_${NET}_AgeMaxSubset_$(date +%F)--$(date +%H)":"$(date +%M).out

# Compute correlations from Network predictions
Rscript --vanilla Scripts/${WDIR}/Compute_NetworkCorrelations.R "sampled_centered" "Outputs/3_GTExDataPrep/Subset_Data" Outputs/${WDIR}/${SUBDIR} "Spearman" &> Logs/${WDIR}/${SUBDIR}/Compute_NetworkCorrelations_${NET}_sampled_centered_Pearson_$(date +%F)--$(date +%H)":"$(date +%M).out
# Rscript --vanilla Scripts/${WDIR}/Compute_NetworkCorrelations.R "sampled_centered" "Outputs/3_GTExDataPrep/Subset_Data/Max_Subset" Outputs/${WDIR}/${SUBDIR}/Age_MaxSubset "Spearman" &> Logs/${WDIR}/${SUBDIR}/Compute_NetworkCorrelations_${NET}_AgeMaxSubset_sampled_centered_Pearson_$(date +%F)--$(date +%H)":"$(date +%M).out

# Compute correlation fold changes
Rscript --vanilla Scripts/${WDIR}/Compute_AgeSlopes.R "sampled_centered" Outputs/${WDIR}/${SUBDIR} &> Logs/${WDIR}/${SUBDIR}/Compute_AgeSlopes_${NET}_AgeSubset_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Compute_AgeSlopes.R "sampled_centered" Outputs/${WDIR}/${SUBDIR}/Age_MaxSubset

# Compute background distribution for correlation fold changes
Rscript --vanilla Scripts/${WDIR}/Compute_AgeSlopeBG.R "sampled_centered" Outputs/${WDIR}/${SUBDIR} &> Logs/${WDIR}/${SUBDIR}/Compute_AgeSlopeBG_${NET}_AgeSubset_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Compute_AgeSlopeBG.R "sampled_centered" Outputs/${WDIR}/${SUBDIR}/Age_MaxSubset

# Heatmap of slopes
Rscript --vanilla Scripts/${WDIR}/Plot_SlopeHeatmap.R

# Volcano plot of slopes
Rscript --vanilla Scripts/${WDIR}/Plot_SlopeVolcano.R

# GO enrichment
Rscript --vanilla Scripts/${WDIR}/Compute_GOEnrichment.R

# Create figure panels
Rscript --vanilla Scripts/${WDIR}/Generate_Plots.R Outputs/${WDIR}/${SUBDIR}/Age_MaxSubset "main" T # Figure 2
Rscript --vanilla Scripts/${WDIR}/Generate_Plots.R Outputs/${WDIR}/${SUBDIR} "all" T # Figure S3
Rscript --vanilla Scripts/${WDIR}/Generate_Plots.R Outputs/5_Predictability/Age/Age_MaxSubset "main" F # Figure S4
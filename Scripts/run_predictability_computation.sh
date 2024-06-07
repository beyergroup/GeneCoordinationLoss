#!/bin/bash

NET=$1
WDIR=$2
echo "Mode = ${NET}"

SUBDIR="WellPredicted_TissueFilters"
mkdir -p Outputs/${WDIR}/${SUBDIR}
mkdir -p Plots/${WDIR}/${SUBDIR}
mkdir -p Logs/${WDIR}/${SUBDIR}

# Randomized network predictions per tissue
Rscript --vanilla Scripts/${WDIR}/Predict_Expression.R ${NET}_randomized "sampled_centered" "Outputs/3_GTExDataPrep/Subset_Data/CrossAge" Outputs/${WDIR}/${SUBDIR}/Randomized &> Logs/${WDIR}/${SUBDIR}/Predict_Expression_${NET}_randomized_$(date +%F)--$(date +%H)":"$(date +%M).out

# Compute correlations from randomized network predictions
Rscript --vanilla Scripts/${WDIR}/Compute_NetworkCorrelations.R "sampled_centered" "Outputs/3_GTExDataPrep/Subset_Data/CrossAge" Outputs/${WDIR}/${SUBDIR}/Randomized "Pearson" &> Logs/${WDIR}/${SUBDIR}/Compute_NetworkCorrelations_${NET}_randomized_sampled_centered_Pearson_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Compute_NetworkCorrelations.R "sampled_centered" "Outputs/3_GTExDataPrep/Subset_Data/CrossAge" Outputs/${WDIR}/${SUBDIR}/Randomized "Spearman" &> Logs/${WDIR}/${SUBDIR}/Compute_NetworkCorrelations_${NET}_randomized_sampled_centered_Spearman_$(date +%F)--$(date +%H)":"$(date +%M).out

# True Network predictions per tissue
Rscript --vanilla Scripts/${WDIR}/Predict_Expression.R ${NET} "sampled_centered" "Outputs/3_GTExDataPrep/Subset_Data/CrossAge" Outputs/${WDIR}/${SUBDIR} &> Logs/${WDIR}/Predict_Expression_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

# Compute correlations from true network predictions
Rscript --vanilla Scripts/${WDIR}/Compute_NetworkCorrelations.R "sampled_centered" "Outputs/3_GTExDataPrep/Subset_Data/CrossAge" Outputs/${WDIR}/${SUBDIR} "Pearson" &> Logs/${WDIR}/${SUBDIR}/Compute_NetworkCorrelations_${NET}_sampled_centered_Pearson_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Compute_NetworkCorrelations.R "sampled_centered" "Outputs/3_GTExDataPrep/Subset_Data/CrossAge" Outputs/${WDIR}/${SUBDIR} "Spearman" &> Logs/${WDIR}/${SUBDIR}/Compute_NetworkCorrelations_${NET}_sampled_centered_Spearman_$(date +%F)--$(date +%H)":"$(date +%M).out


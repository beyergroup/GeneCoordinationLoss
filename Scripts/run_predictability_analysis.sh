#!/bin/bash

NET=$1
WDIR=$2
echo "Mode = ${NET}"

mkdir Outputs/${WDIR}
mkdir Plots/${WDIR}
mkdir Logs/${WDIR}


# Compute well predicted genes based on random distribution

SUBDIR="WellPredicted_TissueFilters"
mkdir Outputs/${WDIR}/${SUBDIR}
mkdir Plots/${WDIR}/${SUBDIR}
mkdir Logs/${WDIR}/${SUBDIR}

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

# Compute correlation thresholds
Rscript --vanilla Scripts/${WDIR}/Compute_WellPredicted.R "sampled_centered_correlations" Outputs/${WDIR}/${SUBDIR} Outputs/${WDIR}/${SUBDIR}/Randomized "Spearman" &> Logs/${WDIR}/${SUBDIR}/Compute_WellPredicted_${NET}_sampled_centered_Spearman_correlations_$(date +%F)--$(date +%H)":"$(date +%M).out

# Characterize well predicted genes



# GTEx tissues and TS cell types



# Age-related changes in predictability

SUBDIR="Age"
mkdir Outputs/${WDIR}/${SUBDIR}
mkdir Plots/${WDIR}/${SUBDIR}
mkdir Logs/${WDIR}/${SUBDIR}

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

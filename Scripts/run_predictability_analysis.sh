
# Compute well predicted genes in a tissue-specific manner

# Randomized Network predictions per tissue
RScript --vanilla PredictExpression.R "stabsel_filtered_largestCC" "GTEx_Networks/Tissue_Networks/Outputs" "Outputs/5_Predictability/WellPredicted_TissueFilters/Randomized"

# Randomized Network predictions per tissue
RScript --vanilla PredictExpression.R "stabsel_filtered_largestCC" "GTEx_Networks/Tissue_Networks/Outputs" "Outputs/5_Predictability/WellPredicted_TissueFilters/Randomized"

# Compute correlations from randomized Network predictions
RScript --vanilla ComputeNetworkCorrelations.R "GTEx_Networks/Tissue_Networks/Outputs" "Outputs/5_Predictability/WellPredicted_TissueFilters/Randomized"

# Compute correlation thresholds
RScript --vanilla Compute_CorrelationThresholds.R "Outputs/5_Predictability/WellPredicted_TissueFilters" "Outputs/5_Predictability/WellPredicted_TissueFilters/Randomized"


# True Network predictions per age group and tissue
RScript --vanilla PredictExpression.R "stabsel_filtered_largestCC" "GTEx_Networks/AgeTissue_Networks/Outputs" "Outputs/5_Predictability/"

# Compute correlations from Network predictions
RScript --vanilla ComputeNetworkCorrelations.R "GTEx_Networks/AgeTissue_Networks/Outputs" "Outputs/5_Predictability"

# Compute correlation fold changes

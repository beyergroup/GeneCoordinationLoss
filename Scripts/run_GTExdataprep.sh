#!/bin/bash

DATA=$1
WDIR=$2
echo "Mode = ${DATA}"

mkdir Outputs/${WDIR}
mkdir Plots/${WDIR}
mkdir Logs/${WDIR}
mkdir Tables/${WDIR}

# Subset per tissue and age group
Rscript --vanilla Scripts/${WDIR}/Subset_Age.R Outputs/GTEx_Preprocessing Outputs/${WDIR}/Subset_Data Plots/${WDIR} &> Logs/${WDIR}/Subset_Age_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Subset_Age_Maximal.R Outputs/GTEx_Preprocessing Outputs/${WDIR}/Subset_Data/Max_Subset &> Logs/${WDIR}/Subset_Age_Maximal_$(date +%F)--$(date +%H)":"$(date +%M).out

# Subset per tissue
Rscript --vanilla Scripts/${WDIR}/Subset_CrossAge.R Outputs/GTEx_Preprocessing Outputs/${WDIR}/Subset_Data Outputs/${WDIR}/Subset_Data/CrossAge &> Logs/${WDIR}/Subset_CrossAge_$(date +%F)--$(date +%H)":"$(date +%M).out

# Center
Rscript --vanilla Scripts/${WDIR}/Center_NormalizedData.R Outputs/${WDIR}/Subset_Data &> Logs/${WDIR}/Center_Subset_Age_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Center_NormalizedData.R Outputs/${WDIR}/Subset_Data/Max_Subset &> Logs/${WDIR}/Center_Subset_Age_Maximal_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Center_NormalizedData.R Outputs/${WDIR}/Subset_Data/CrossAge &> Logs/${WDIR}/Center_Subset_CrossAge_$(date +%F)--$(date +%H)":"$(date +%M).out

# Count samples per subset (Supplementary Tables 1 and 2)
Rscript --vanilla Scripts/${WDIR}/Count_Samples.R Outputs/${WDIR}/Subset_Data Tables/${WDIR}/TableS2.xlsx &> Logs/${WDIR}/Count_Subset_Age_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Count_Samples.R Outputs/${WDIR}/Subset_Data/Max_Subset Tables/${WDIR}/TableS1.xlsx &> Logs/${WDIR}/Count_Subset_Age_Maximal_$(date +%F)--$(date +%H)":"$(date +%M).out
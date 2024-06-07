#!/bin/bash

DATA=$1
WDIR=$2
echo "Mode = ${DATA}"

mkdir Outputs/${WDIR}
mkdir Plots/${WDIR}
mkdir Logs/${WDIR}

# Load and arrange data
Rscript --vanilla Scripts/${WDIR}/Arrange_Data.R "Data/GTEx" Outputs/${WDIR} &> Logs/${WDIR}/Arrange_Data_$(date +%F)--$(date +%H)":"$(date +%M).out

# Normalize with DESeq2
Rscript --vanilla Scripts/${WDIR}/Normalize_DESeq2.R Outputs/${WDIR} &> Logs/${WDIR}/Normalize_DESeq2_$(date +%F)--$(date +%H)":"$(date +%M).out

# Filter samples
Rscript --vanilla Scripts/${WDIR}/Filter_Samples.R Outputs/${WDIR} &> Logs/${WDIR}/Filter_Samples_$(date +%F)--$(date +%H)":"$(date +%M).out

# Filter genes
Rscript --vanilla Scripts/${WDIR}/Filter_Genes.R Outputs/${WDIR} Plots/${WDIR} &> Logs/${WDIR}/Filter_Genes_$(date +%F)--$(date +%H)":"$(date +%M).out

# Perform batch correction
Rscript --vanilla Scripts/${WDIR}/Correct_BatchEffects.R Outputs/${WDIR} &> Logs/${WDIR}/Correct_BatchEffects_$(date +%F)--$(date +%H)":"$(date +%M).out

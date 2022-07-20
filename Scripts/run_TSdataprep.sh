#!/bin/bash

WDIR=$1

mkdir Outputs/${WDIR}
mkdir Plots/${WDIR}
mkdir Logs/${WDIR}


# Unzip external data
unzip External/Tabula_Sapiens/*.zip

# Arrange data
Rscript --vanilla Scripts/${WDIR}/Arrange_Data.R "External/Tabula_Sapiens/" ${WDIR} &> Logs/${WDIR}/Arrange_Data_$(date +%F)--$(date +%H)":"$(date +%M).out

# Remove low quality cells
Rscript --vanilla Scripts/${WDIR}/Filter_Cells.R ${WDIR} "10X" &> Logs/${WDIR}/Filter_Cells_10X_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Filter_Cells.R ${WDIR} "smartseq2" &> Logs/${WDIR}/Filter_Cells_smartseq2_$(date +%F)--$(date +%H)":"$(date +%M).out

# Remove unreliable genes
Rscript --vanilla Scripts/${WDIR}/Filter_Genes.R ${WDIR} "10X" &> Logs/${WDIR}/Filter_Genes_10X_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Filter_Genes.R ${WDIR} "smartseq2" &> Logs/${WDIR}/Filter_Genes_smartseq2_$(date +%F)--$(date +%H)":"$(date +%M).out

# Normalize


# Subset
Rscript --vanilla Scripts/${WDIR}/Subset_SC.R ${WDIR} "10X" &> Logs/${WDIR}/Subset_SC_10X_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Subset_SC.R ${WDIR} "smartseq2" &> Logs/${WDIR}/Subset_SC_smartseq2_$(date +%F)--$(date +%H)":"$(date +%M).out

# Center
Rscript --vanilla Scripts/${WDIR}/Center_NormalizedData.R ${WDIR} &> Logs/${WDIR}/Center_NormalizedData_$(date +%F)--$(date +%H)":"$(date +%M).out

# Quantile transform
Rscript --vanilla Scripts/${WDIR}/Transform_Quantile.R ${WDIR} "10X" &> Logs/${WDIR}/Transform_Quantile_10X_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Transform_Quantile.R ${WDIR} "smartseq2" &> Logs/${WDIR}/Transform_Quantile_smartseq2_$(date +%F)--$(date +%H)":"$(date +%M).out

# Create cross-tissue with cell types that passed all steps



# Create pseudobulk per cell type and individual
Rscript --vanilla Scripts/${WDIR}/Create_Pseudobulk.R ${WDIR} "10X" &> Logs/${WDIR}/Create_Pseudobulk_10X_$(date +%F)--$(date +%H)":"$(date +%M).out
Rscript --vanilla Scripts/${WDIR}/Create_Pseudobulk.R ${WDIR} "smartseq2" Logs/${WDIR}/Create_Pseudobulk_smartseq2_$(date +%F)--$(date +%H)":"$(date +%M).out


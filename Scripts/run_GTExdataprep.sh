#!/bin/bash

DATA=$1
WDIR=$2
echo "Mode = ${DATA}"

mkdir Outputs/${WDIR}
mkdir Plots/${WDIR}
mkdir Logs/${WDIR}


# Compute GTEx tissue LFCs

# Compute GTEx tissue variances

# Create pseudobulk for Tabula Sapiens
Rscript --vanilla Create_Pseudobulk.R $


# Single cell data
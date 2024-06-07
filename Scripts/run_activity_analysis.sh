#!/bin/bash

NET=$1
WDIR=$2
echo "Mode = ${NET}"

mkdir Outputs/${WDIR}
mkdir Plots/${WDIR}
mkdir Logs/${WDIR}


# Compute GTEx tissue LFCs

# Compute GTEx tissue variances

# Smooth GTEx LFCs through all modules (GTEx)
Rscript --vanilla Scripts/${WDIR}/Smooth_LFCs.R ${NET} ${WDIR} "rownormalized" "" &> Logs/${WDIR}/Smooth_LFCs_RWR_${NET}_rownormalized_$(date +%F)--$(date +%H)":"$(date +%M).out

# Average smoothed values across gene modules (GTEx)
Scripts/${WDIR}/Average_SmoothedLFCs.R

# Plot obtained relative activities (GTEx)
Scripts/${WDIR}/Plot_RelativeActivities.R


# Average variance, within-cluster correlation


# Single cell data
#!/bin/bash

NET=$1
WDIR=$2
echo "Mode = ${NET}"

mkdir Outputs/${WDIR}
mkdir Plots/${WDIR}
mkdir Logs/${WDIR}

# Collect information on genomic position of genes
Rscript --vanilla Scripts/${WDIR}/Gather_GeneInfo.R &> Logs/${WDIR}/Gather_GeneInfo_$(date +%F)--$(date +%H)":"$(date +%M).out

# Compute moving average of age slopes
Rscript --vanilla Scripts/${WDIR}/Compute_MovingAverage.R &> Logs/${WDIR}/Compute_MovingAverage_$(date +%F)--$(date +%H)":"$(date +%M).out


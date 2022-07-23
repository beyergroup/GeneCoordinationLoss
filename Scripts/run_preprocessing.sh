#!/bin/bash

NET=$1
WDIR=$2
echo "Mode = ${NET}"

mkdir Outputs/${WDIR}
mkdir Plots/${WDIR}
mkdir Logs/${WDIR}

# Look at in- vs out-going edges
echo "Plotting in- and out-going edge comparisons"
Rscript --vanilla Scripts/${WDIR}/Plot_InOutEdges.R ${NET} ${WDIR} &> Logs/${WDIR}/Plot_InOutEdges_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

# Remove unbalanced edges (should be directed but are not):
echo "Removing unbalanced edge pairs"
Rscript --vanilla Scripts/${WDIR}/Remove_UnbalancedEdgePairs.R ${NET} ${WDIR} &> Logs/${WDIR}/Remove_UnbalancedEdgePairs_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

# Remove predictors located in same chromosome arm:
echo "Removing cis-acting predictors"
Rscript --vanilla Scripts/${WDIR}/Remove_CisPredictors.R ${NET}_filtered ${WDIR} &> Logs/${WDIR}/Remove_CisPredictors_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

# Reduce to largest connected component:
echo "Reducing network to largest connected component"
Rscript --vanilla Scripts/${WDIR}/Reduce_ToLargestCC.R ${NET}_filtered_trans ${WDIR} &> Logs/${WDIR}/Reduce_ToLargestCC_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

#!/bin/bash

NET=$1
WDIR=$2
echo "Mode = ${NET}"

mkdir -p Outputs/${WDIR}
mkdir -p Plots/${WDIR}
mkdir -p Logs/${WDIR}
mkdir -p Tables/${WDIR}

# Look at in- vs out-going edges
echo "Plotting in- and out-going edge comparisons"
Rscript --vanilla Scripts/${WDIR}/Plot_InOutEdges.R ${NET} ${WDIR} &> Logs/${WDIR}/Plot_InOutEdges_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

# Remove unbalanced edges (should be directed but are not):
echo "Removing unbalanced edge pairs"
Rscript --vanilla Scripts/${WDIR}/Remove_UnbalancedEdgePairs.R ${NET} ${WDIR} &> Logs/${WDIR}/Remove_UnbalancedEdgePairs_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

# Remove predictors located in same chromosome arm:
echo "Removing cis-acting predictors"
Rscript --vanilla Scripts/${WDIR}/Remove_CisPredictors.R ${NET}_filtered ${WDIR} &> Logs/${WDIR}/Remove_CisPredictors_${NET}_filtered_$(date +%F)--$(date +%H)":"$(date +%M).out

# Reduce to largest connected component:
echo "Reducing network to largest connected component"
Rscript --vanilla Scripts/${WDIR}/Reduce_ToLargestCC.R ${NET}_filtered_trans ${WDIR} &> Logs/${WDIR}/Reduce_ToLargestCC_${NET}_filtered_trans_$(date +%F)--$(date +%H)":"$(date +%M).out

# Characterize network:
echo "Characterizing network"
Rscript --vanilla Scripts/${WDIR}/Characterize.R ${NET}_filtered_trans_largestCC ${WDIR} &> Logs/${WDIR}/Characterize_${NET}_filtered_trans_$(date +%F)--$(date +%H)":"$(date +%M).out

# Randomize network:
echo "Randomizing network edges for negative control"
Rscript --vanilla Scripts/${WDIR}/Randomize_Network.R ${NET}_filtered_trans_largestCC ${WDIR} &> Logs/${WDIR}/Randomize_Network_${NET}_filtered_trans_largestCC_$(date +%F)--$(date +%H)":"$(date +%M).out

#!/bin/bash

NET=$1
WDIR=$2
echo "Mode = ${NET}"

mkdir Outputs/${WDIR}
mkdir Plots/${WDIR}
mkdir Logs/${WDIR}

# Look at in- vs out-going edges
Rscript --vanilla Scripts/${WDIR}/Plot_InOutEdges.R ${NET} ${WDIR} &> Logs/${WDIR}/Plot_InOutEdges_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

# Remove unbalanced edges (should be directed but are not):
Rscript --vanilla Scripts/${WDIR}/Remove_UnbalancedEdgePairs.R ${NET} ${WDIR} &> Logs/${WDIR}/Remove_UnbalancedEdgePairs_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

# Reduce to largest connected component:
Rscript --vanilla Scripts/${WDIR}/Reduce_ToLargestCC.R ${NET}_filtered ${WDIR} &> Logs/${WDIR}/Reduce_ToLargestCC_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out


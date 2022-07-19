#!/bin/bash

NET=$1
WDIR=$2
echo "Mode = ${NET}"

mkdir Outputs/${WDIR}
mkdir Plots/${WDIR}
mkdir Logs/${WDIR}


# Compute centrality metrics
Rscript --vanilla Scripts/${WDIR}/Compute_Centrality.R ${NET} ${WDIR} &> Logs/${WDIR}/Compute_Centrality_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

# Plot centrality metrics without comparison to poorly predicted genes
Rscript --vanilla Scripts/${WDIR}/Plot_Centrality.R ${NET} ${WDIR} "FALSE" &> Logs/${WDIR}/Plot_Centrality_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

# GO term enrichment on centrality measures (compute & plot)
Rscript --vanilla Scripts/${WDIR}/Run_GOEnrichment.R ${NET} ${WDIR} "FALSE" &> Logs/${WDIR}/Run_GOEnrichment_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out
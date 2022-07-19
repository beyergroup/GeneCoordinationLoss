#!/bin/bash

NET=$1
WDIR=$2
echo "Mode = ${NET}"

mkdir Outputs/${WDIR}
mkdir Plots/${WDIR}
mkdir Logs/${WDIR}


# Get network representation (adjacency or Laplacian matrices)
Rscript --vanilla Scripts/${WDIR}/Compute_NetworkRepresentation.R ${NET} ${WDIR} &> Logs/${WDIR}/Compute_NetworkRepresentation_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

# Decompose into eigenvectors
Rscript --vanilla Scripts/${WDIR}/Decompose_Network.R ${NET} ${WDIR} "undirected" &> Logs/${WDIR}/Decompose_Network_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

# Hierarchical clustering of eigenvectors
Rscript --vanilla Scripts/${WDIR}/Cluster_Eigenvectors.R ${NET} ${WDIR} "rownormalized" &> Logs/${WDIR}/Cluster_Eigenvectors_${NET}_rownormalized_$(date +%F)--$(date +%H)":"$(date +%M).out

# Compute modularity to find optimal clustering decision
Rscript --vanilla Scripts/${WDIR}/Compute_Modularity.R ${NET} ${WDIR} "rownormalized" &> Logs/${WDIR}/Compute_Modularity_${NET}_rownormalized_$(date +%F)--$(date +%H)":"$(date +%M).out

# Plot modularity heatmaps
Rscript --vanilla Scripts/${WDIR}/Plot_Modularity.R ${NET} ${WDIR} "rownormalized" &> Logs/${WDIR}/Plot_Modularity_${NET}_rownormalized_$(date +%F)--$(date +%H)":"$(date +%M).out

# Compute modules based on optimal clustering decision
Rscript --vanilla Scripts/${WDIR}/Call_Modules.R ${NET} ${WDIR} "rownormalized" &> Logs/${WDIR}/Call_Modules_${NET}_rownormalized_$(date +%F)--$(date +%H)":"$(date +%M).out

# Characterize the computed modules
Rscript --vanilla Scripts/${WDIR}/Plot_ModuleFeatures.R ${NET} ${WDIR} "rownormalized" &> Logs/${WDIR}/Plot_ModuleFeatures_${NET}_rownormalized_$(date +%F)--$(date +%H)":"$(date +%M).out

# GO enrichment in computed modules
Rscript --vanilla Scripts/${WDIR}/Run_GOEnrichment.R ${NET} ${WDIR} "rownormalized" &> Logs/${WDIR}/Run_GOEnrichment_${NET}_rownormalized_$(date +%F)--$(date +%H)":"$(date +%M).out

# Plot GO enrichment of selected modules
Rscript --vanilla Scripts/${WDIR}/Plot_GOEnrichment.R ${NET} ${WDIR} "sum_rownormalized" "1" &> Logs/${WDIR}/Plot_GOEnrichment_${NET}_sum_rownormalized_1_$(date +%F)--$(date +%H)":"$(date +%M).out

# At this point the decision is made for sum of weights, so we further cluster the big module obtained with that analysis

#!/bin/bash

#SBATCH --account=adesous1
#SBATCH --partition=single
#SBATCH --ntasks=1
#SBATCH --exclusive

echo $threshold
Rscript --vanilla compute_totaleffects_network.R "/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel/Networks/Adjacency_directed_weightsident_rownormalized_largest_cc.rds" $threshold "Outputs/Human_Network/stabsel/Networks/TotalEffects/"
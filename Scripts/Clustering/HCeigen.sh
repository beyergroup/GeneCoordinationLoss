#!/bin/bash

#SBATCH --account=adesous1
#SBATCH --error=../../Logs/HCeigen_unweightedAdjacency.err
#SBATCH --output=../../Logs/HCeigen_unweightedAdjacency.out
#SBATCH --job-name=HCeigen
#SBATCH --partition=single
#SBATCH --nodelist=beyer-n02
#SBATCH --ntasks=1

Rscript --vanilla hc_eigenvectors.R "stabsel" "network_largest_cc" "none" "Adjacency"
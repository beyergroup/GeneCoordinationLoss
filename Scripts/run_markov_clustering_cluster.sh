#!/bin/bash

#SBATCH --account=adesous1
#SBATCH --error=../Logs/hbmMCL_cluster_%j.err
#SBATCH --output=../Logs/hbmMCL_cluster_%j.out
#SBATCH --job-name=hbmMarkovCL
#SBATCH --partition=single
#SBATCH --nodelist=beyer-n03
#SBATCH --ntasks=1

module unload R-3.5.1
module load R-3.4.4

Rscript --vanilla markov_clustering_cluster.R stabsel
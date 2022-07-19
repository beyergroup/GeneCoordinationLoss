#!/bin/bash
#SBATCH --account=adesous1
#SBATCH --partition=single
#SBATCH -c 1

module unload R-3.5.1
module load R-3.4.4

echo $tissue
Rscript --vanilla ZINBWaVe_normalization.R "../../Outputs/Tabula_Sapiens/GeneFiltered/all_detected_genes" "../../Outputs/Tabula_Sapiens/ZBNormalized/all_detected_genes" TRUE $tissue
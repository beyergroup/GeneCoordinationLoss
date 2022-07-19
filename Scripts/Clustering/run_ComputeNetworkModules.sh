#!/bin/bash

MAINDIR="/data/public/adesous1/GeneCorrelation/"
INDIR=${MAINDIR}"Outputs/Human_Network/stabsel/Networks"
files=$(ls $INDIR *undirected*)
net="stabsel"
analysis="Modules/Eigen"

for file in $files
do
  echo $file
  Rscript --vanilla ComputeNetworkModules.R $file ${MAINDIR}"Outputs/Human_Network/"${net}"/Networks/" ${MAINDIR}"Outputs/Human_Network/"${net}"/"${analysis}"/" ${MAINDIR}"Plots/Human_Network/"${net}"/"${analysis}"/" >> ${MAINDIR}"Logs/ComputeNetworkModules_"$file".out"
done
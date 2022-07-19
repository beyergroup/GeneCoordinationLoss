#!/bin/bash -l

nodes=(beyer-n01 beyer-n02 beyer-n03)
tissues=(Lung Pancreas Vasculature)

for i in 0 1 2;
do
echo ${tissues[$i]}
echo ${nodes[$i]}
export tissue=${tissues[$i]}
sbatch -n1 -N1 -w ${nodes[$i]} --error=ZINBWaVe_cluster_${i}.log --job-name=ZB$i  ZINBWaVe_normalization_cluster.sh &
done

wait


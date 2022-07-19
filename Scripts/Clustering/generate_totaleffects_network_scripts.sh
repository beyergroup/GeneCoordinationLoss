#!/bin/bash -l

thre=(1 0.1 0.01)

for i in 0 1 2
do
  echo ${thre[$i]}
  export threshold=${thre[$i]}
  sbatch -n1 -N1 --error=../../Logs/MP_${threshold}.log --output=../../Logs/MP_${threshold}.out --job-name=MP_${threshold}  run_totaleffects_network.sh &
done
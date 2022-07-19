#!/bin/bash -l

thre=(1 0.1)

for i in 1
do
  echo ${thre[$i]}
  export threshold=${thre[$i]}
  sbatch -n1 -N1 --error=../../Logs/MP_${threshold}_test.log --output=../../Logs/MP_${threshold}_test.out --job-name=MP_${threshold}_test  run_tmp_totaleffects_test.sh &
done
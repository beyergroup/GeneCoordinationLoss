#!/bin/bash

INDIR=$1
NET=$2
OUTDIR=$3

Rscript Scripts/Get_Network.R ${INDIR}"/network_funclib.r" ${INDIR}"/${NET}.Rout" ${OUTDIR} &> Logs/Get_Network_${NET}_$(date +%F)--$(date +%H)":"$(date +%M).out

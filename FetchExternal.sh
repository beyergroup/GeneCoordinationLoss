#!/bin/bash

# Grab network handling library
cd /data/public/xwu2/Rscripts/network_funclib.r External/network_funclib.r

# Grab stabsel network
cp /data/public/xwu2/networks/infer_CCLERSEQ_01/Network_Whole_CCLERSEQ_01.Rout External/stabsel.Rout

# Grab stabsel pcclasso network
cp /data/public/xwu2/networks/infer_CCLERSEQ_01/Network_Whole_CCLERSEQ_01.pcc_lasso_regress.regress.Rout External/stabsel_pcclasso.Rout

# Grab randomized network
cp /data/public/xwu2/networks/infer_CCLERSEQ_01/Network_Whole_CCLERSEQ_01.shuffle.per100.Rout External/stabsel_randomized.Rout

# Prediction of expression based on the networks

args <- commandArgs(trailingOnly = T)

NET_FILE = "Data/Networks/Human/stabsel_network_Hs_filtered.rds"
EXPR_DIR = "Outputs/Tabula_Sapiens/SampledSCTNorm/"
OUTDIR = "Outputs/Human_Network/stabsel/Predictability/Tabula_Sapiens/"

setwd("../../")
dir.create(OUTDIR)

source("Scripts/functions.R")

net <- as.matrix(readRDS(NET_FILE))
data_files <- list.files(EXPR_DIR, pattern = "centered_data.rds")

for(file in data_files){
  
  centered_data <- readRDS(paste0(EXPR_DIR,file))
  
  pred <- PredictNet(net, as.matrix(centered_data), maxiter = 1)
  saveRDS(pred, paste0(OUTDIR,
                       gsub("_data","_net_predictions",file)))
  
  rm(centered_data,pred); gc()
}

rm(net,data_files); gc()

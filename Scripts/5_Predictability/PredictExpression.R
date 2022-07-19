# Predict target gene expression

args = commandArgs(trailingOnly = T)
NET = args[1]
EXPR_DIR = args[2]
OUTDIR = args[3]

source("Scripts/functions.R")

net <- as.matrix(ReadRDS(paste0("Outputs/0_Preprocessing/",NET,"_network_Hs.rds")))
data_files <- list.files(EXPR_DIR, pattern = "sampled_centered_data.rds")

for(file in data_files){
  
  centered_data <- readRDS(paste0(EXPR_DIR,"/",file))
  centered_data <- DataENSGToSymbol(centered_data)
  
  pred <- PredictNet(net, as.matrix(centered_data), maxiter = 1)
  saveRDS(pred, paste0(OUTDIR,"/",
                       gsub("_data","_net_predictions",file)))
  
  rm(centered_data,pred); gc()
}

rm(net,data_files); gc()

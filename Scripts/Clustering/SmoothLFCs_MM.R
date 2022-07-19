# Compute Module Activity (considering global topology, iterative multiplication of network without restart)

.libPaths("Resources/Rlibs/R-4.0.3/")
library(devtools)
library(igraph)
source("Scripts/functions.R")

NET_DIR = "/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel/Networks/"
EXPR_DIR = "/data/public/adesous1/GeneCorrelation/Outputs/GTEx/logFC/"
OUT_DIR = "/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel/Modules/Module_Activity/SmoothedLFCs/"

dir.create(OUT_DIR)

net_files <- list.files(NET_DIR)
net_files <- net_files[grep("_rownormalized",net_files)] # only normalized networks
net_files <- net_files[grep("largest_cc",net_files)] # ignore full networks (largest cc is disconnected anyways!)

expr_files <- list.files(EXPR_DIR, pattern = "sampled_meanLFC.rds")


for(net_file in net_files){
  
  network <- readRDS(paste0(NET_DIR,net_file))
  
  smoothed <- list()
  
  for(file in expr_files){
    tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
    LFC <- readRDS(paste0(EXPR_DIR,file))
    smoothed[[tissue]] <- IterativeNetMultiplication(vector = LFC,
                                                     network = network)
  }
  
  smoothed <- do.call(cbind, smoothed)
  rm(file,tissue,LFC); gc()
  
  # save smoothed LFCs
  saveRDS(smoothed, paste0(OUT_DIR,gsub(".rds","_NetMM50.rds",net_file)))
  
  rm(smoothed,network); gc()
}


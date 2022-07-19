# Compute Module Activity (considering global topology, RWR)

args <- commandArgs(trailingOnly = T)
NET = args[1]
WDIR = args[2]
PRED_DIR = args[3]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(devtools)
library(igraph)
install_github("beyergroup/BioNetSmooth")
library(BioNetSmooth)

source("Scripts/functions.R")

NET_DIR = "/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel/Networks/"
PRED_DIR = "/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel/Predictability/AgeTissue/"
OUT_DIR = "/data/public/adesous1/GeneCorrelation/Outputs/5_Predictability/SmoothedLFCs/"

dir.create(OUT_DIR)

net_file <- readRDS("Outputs/0_Preprocessing/stabsel_filtered_largestCC_network_Hs.rds")
predictability_files <- list.files(paste0("Outputs/Human_Network/stabsel/Predictability/AgeTissue"),
                                   pattern = "_sampled_ageDP.rds", full.names = T)


network <- readRDS(paste0(NET_DIR,net_file))

# read in LFCs from all tissues at once
LFC <- sapply(predictability_files, readRDS)
colnames(LFC) <- sapply(expr_files, function(x) strsplit(x,"_")[[1]][1])

  for(alpha in c(0.2,0.5,0.8)){
    
    message("Smoothing for alpha = ",alpha)
    # smooth
    smoothed <- SmoothNetworkRWR(net = network, expr = LFC, alpha)
    
    # save smoothed LFCs
    saveRDS(smoothed, paste0(OUT_DIR, gsub(".rds",paste0("_RWR",alpha,".rds"),
                                           net_file)))
    rm(smoothed); gc()
  }
  rm(LFC,network); gc()
}


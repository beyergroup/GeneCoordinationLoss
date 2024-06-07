# Compute Module Activity (considering global topology, RWR)

args <- commandArgs(trailingOnly = T)
NET = "stabsel_filtered_trans_largestCC"
OUT_DIR = "Outputs/5_Predictability/Age/Age_MaxSubset/SmoothedSlopes"
EXPR_DIR = "Outputs/5_Predictability/Age/Age_MaxSubset"

.libPaths("Resources/Rlibs/R-4.0.3/")
library(devtools)
library(igraph)
install_github("beyergroup/BioNetSmooth")
library(BioNetSmooth)

source("Scripts/functions.R")

NET_DIR = "/data/public/adesous1/GeneCorrelation/Outputs/0_Preprocessing/Adjacency"
network <- ReadRDS(paste0(NET_DIR,"/",NET,"_Adjacency_undirected_weightssum_",
                          "rownormalized.rds"))
network <- MatrixToSquare(network)


dir.create(OUT_DIR)

expr_files <- list.files(EXPR_DIR, pattern = "ageslope_well_predicted.rds")


for(file in expr_files){
  
  slopes <- ReadRDS(paste0(EXPR_DIR,"/",file))
  
  for(alpha in c(0.2,0.5,0.8)){
    
    message("Smoothing for alpha = ",alpha)
    
    # Map values to network
    smoothed <- SmoothNetworkRWR(net = network,
                                 expr = slopes[,c("Slope","pval")],
                                 alpha)
    
    # save smoothed LFCs
    saveRDS(smoothed, paste0(OUT_DIR,"/",gsub(".rds",
                                              paste0("_RWR_",alpha,".rds"),
                                              file)))
    rm(smoothed); gc()
  }
  rm(slopes); gc()
}
rm(LFC,network); gc()


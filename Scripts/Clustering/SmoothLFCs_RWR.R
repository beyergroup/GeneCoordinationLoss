# Compute Module Activity (considering global topology, RWR)
args <- commandArgs(trailingOnly = T)

.libPaths("Resources/Rlibs/R-4.0.3/")
library(devtools)
library(igraph)
install_github("beyergroup/BioNetSmooth")
library(BioNetSmooth)

source("Scripts/functions.R")

NET_DIR = "/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel/Networks/"
EXPR_DIR = "/data/public/adesous1/GeneCorrelation/Outputs/GTEx/logFC/"
OUT_DIR = "/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel/Modules/Module_Activity/SmoothedLFCs/"

TS = T

dir.create(OUT_DIR)
if(TS){
  OUT_DIR <- paste0(OUT_DIR,"Tabula_Sapiens")
  dir.create(OUT_DIR)
}

if(TS){
  net_files <- "Adjacency_directed_weightsident_rownormalized_largest_cc.rds"
  expr_files <- "Outputs/Tabula_Sapiens/Pseudobulk/Mean/pseudobulk.rds"
} else{
  net_files <- list.files(NET_DIR)
  net_files <- net_files[grep("_undirected",net_files)]
  net_files <- net_files[grep("largest_cc",net_files)] # ignore full networks (largest cc is disconnected anyways!)
  net_files <- net_files[-grep("weightsident|signedmaxabs",net_files)] # RWR cannot support negative weights
  
  expr_files <- list.files(EXPR_DIR, pattern = "sampled_meanLFC.rds")
}

for(net_file in net_files){
  
  network <- readRDS(paste0(NET_DIR,net_file))
  
  if(TS){
    LFC <- readRDS(expr_files)
  } else{
    # read in LFCs from all tissues at once
    LFC <- sapply(paste0(EXPR_DIR,expr_files), readRDS)
    colnames(LFC) <- sapply(expr_files, function(x) strsplit(x,"_")[[1]][1])
  }
  
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


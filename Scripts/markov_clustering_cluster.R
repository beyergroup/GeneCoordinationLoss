# Markov Clustering

setwd("../")

args = commandArgs(trailingOnly=TRUE)
net = args[1]
netfile = paste0("Outputs/Human_Network/",net,"/undirected_adj_mat.rds")

library(Matrix)
# library(expm, lib.loc = "Resources/Rlibs/R-3.4.4/")
# source("Resources/Rlibs/MCL/R/mcl.R")
library(hbm, lib.loc = "Resources/Rlibs/R-3.4.4/")

adj_undirected <- readRDS(netfile)
adj_undirected <- as.matrix(adj_undirected)
adj_undirected <- abs(adj_undirected)

for(infl in c(1.5,2,3,4,5)){
  # out <- mcl(adj_undirected, inflation = infl, addLoops = F, max.iter = 50)
  out <- mcl(adj_undirected, infl = infl)
  saveRDS(out, paste0("Outputs/Human_Network/",net,"/Topology/Modules/hbmMCL_",infl,"_out.rds"))
}


# ClusterONE modules

source("/cellnet/GeneCorrelation/Human_Network_characterization/functions.R")
.libPaths(c("/data/public/adesous1/GeneCorrelation/Resources/Rlibs/R-4.0.3/",.libPaths()))
library(igraph)
library(topGO)

setwd("../")

net_graph <- readRDS("Outputs/Human_Network/network_igraph.rds")

clusters <- read_csv("Outputs/Human_Network/Topology/ClusterONE/0.3/clusters.csv")

for(i in 1:nrow(clusters)){
  genes <- strsplit(gsub("\"", "", clusters$Members[i]), " ")[[1]]
  GetGOEnrich(genes, names(V(net_graph)), "BP", enrich_cutoff = 0, pval_cutoff = 0, algorithm = "classic")
  GetGOEnrich(genes, names(V(net_graph)), "MF", enrich_cutoff = 0, pval_cutoff = 0, algorithm = "classic")
}

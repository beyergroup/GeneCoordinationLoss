source("/cellnet/GeneCorrelation/Human_Network_characterization/functions.R")
.libPaths(c("/data/public/adesous1/GeneCorrelation/Resources/Rlibs/R-4.0.3/",.libPaths()))
library(igraph)
library(topGO)

setwd("../")

net_graph <- readRDS("Outputs/Human_Network/Topology/network_igraph.rds")
# remove directionality to find communities
net_graph_undirected <- as.undirected(net_graph, mode = "collapse", edge.attr.comb = "mean")

edge_list <- as_edgelist(net_graph_undirected)
write.table(edge_list, "Outputs/Human_Network/undirected_edge_list.txt", sep = "\t",
            row.names = F, col.names = F)

edge_list <- get.data.frame(net_graph_undirected)
edge_list$weight <- abs(edge_list$weight)
edge_list$weight <- edge_list$weight/(max(edge_list$weight))
write.table(edge_list, "Outputs/Human_Network/undirected_abs_scale_weights_edge_list.txt", sep = "\t",
            row.names = F, col.names = F)

# Get adjacency and Laplacian matrices

args <- commandArgs(trailingOnly = T)
NET = args[1]
WDIR = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(igraph)
source("Scripts/functions.R")


params <- expand.grid("Directed" = c(TRUE,FALSE),
                      "Weights" = c("ident","abs","none","maxabs",
                                    "signedmaxabs","sum"),
                      "Normalized" = c(TRUE,FALSE),
                      stringsAsFactors = F)
params <- subset(params, (Directed & (Weights %in% c("ident","abs","none"))) |
                   (!Directed & (Weights %in% c("maxabs","signedmaxabs",
                                                "sum","none"))))

# Read in original network
net_graph <- ReadRDS(paste0("Outputs/0_Preprocessing/",NET,"_network_Hs.rds"))
net_graph <- igraph::graph_from_adjacency_matrix(net_graph, mode = "directed",
                                                 weighted = T)

for(i in 1:nrow(params)){
  CreateNetworkMatrix(param = params[i,],
                      net_graph = net_graph,
                      path = paste0("Outputs/0_Preprocessing/Adjacency/",NET,"_"))
}


.libPaths("Resources/Rlibs/R-4.0.3/")
library(igraph)
source("Scripts/functions.R")


# Adjacency matrix - directed: original, abs or no weights
params <- data.frame("Full" = T, "Directed" = T, "Weights" = c("ident","abs","none"),
                     "Normalized" = F)

# Adjacency matrix - directed, largest connected component: original, abs or no weights
params <- rbind.data.frame(params,
                           data.frame("Full" = F, "Directed" = T,
                                      "Weights" = c("ident","abs","none"),
                                      "Normalized" = F))

# Adjacency matrix - undirected, largest connected component: max abs, max signed abs or no weights
params <- rbind.data.frame(params,
                           data.frame("Full" = F, "Directed" = F,
                                      "Weights" = c("maxabs","signedmaxabs","none"),
                                      "Normalized" = F))

# Row-normalized versions
tmp <- params
tmp[,"Normalized"] <- T
params <- rbind.data.frame(params, tmp)
rm(tmp)


# Read in original network
net_graph <- readRDS("Outputs/Human_Network/stabsel/network_igraph.rds")

for(i in 1:nrow(params)){
  CreateNetworkMatrix(param = params[i,],
                      net_graph = net_graph,
                      path = "Outputs/Human_Network/stabsel/Networks/")
}

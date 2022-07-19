# Components of the network (no PCC)

source("functions.R")

# No removal of residual edges
net <- readRDS("../Data/Networks/Human/stabsel_network_Hs.rds")
net <- as.matrix(net)
net_sq <- MatrixToSquare(net)

library(igraph, lib.loc = "../Resources/Rlibs/R-4.0.3/")
g <- graph_from_adjacency_matrix(net_sq, mode = "directed", weighted = T)
c <- components(g) # 22050 genes in largest connected component

# Undirected with both edges taking highest absolute value
g <- as.undirected(g, mode = "collapse", 
                   edge.attr.comb = list(weight = function(x) x[which.max(abs(x))]))
c <- components(g) # 22050 genes in largest connected component too

# Removal of residual edges
net <- readRDS("../Data/Networks/Human/stabsel_network_Hs_filtered.rds")

gg <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
cc <- components(gg) # 22050 genes in largest connected component too.

rm(net, net_sq, g, c, gg, cc); gc()

# Removal of residual edges ofc makes no difference bc we are only removing an
# edge when a bigger one exists in the opposite direction...


# Components of trimmed network

net <- readRDS("../Data/Networks/Human/stabsel_pcclasso_network_Hs.rds")
net <- as.matrix(net)

g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
c <- components(g) # 14607 genes in largest connect component

rm(net, g, c); gc()
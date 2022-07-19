# Reduce to largest connected component

args <- commandArgs(trailingOnly = T)
NET = args[1]
WDIR = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")
library(igraph)

# Load network ----------------------------------------------------------------
net <- ReadRDS(paste0("Outputs/",WDIR,"/",NET,"_network_Hs.rds"))
net <- as.matrix(net)
net <- MatrixToSquare(net)

# Convert to igraph -----------------------------------------------------------
g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
c <- components(g)

message(max(c$csize)," genes in largest connected component")

# Subset network --------------------------------------------------------------
vids <- V(g)[c$membership == which.max(c$csize)]
sub_g <- induced_subgraph(g, vids = vids)
rm(g,c,net); gc()

# Save output -----------------------------------------------------------------
sub_net <- as_adjacency_matrix(sub_g, type = "both", attr = "weight")
WriteRDS(sub_net, paste0("Outputs/",WDIR,"/",NET,"_largestCC_network_Hs.rds"))

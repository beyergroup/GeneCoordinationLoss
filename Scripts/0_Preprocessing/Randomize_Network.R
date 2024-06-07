# Randomize network

args = commandArgs(trailingOnly = T)
NET = args[1]
WDIR = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")

net <- as.matrix(ReadRDS(paste0("Outputs/",WDIR,"/",NET,"_network_Hs.rds")))

# retrieve edges (predictors on 1st col and targets on 2nd col)
edges <- which(net != 0, arr.ind = T)
# add weight information
edges <- cbind(edges,
               "Weight" = apply(edges, 1, function(i) net[i[1],i[2]]))

# Randomly shuffle predictors, while preserving targets and weights
set.seed(1)
shuffled_edges <- edges
shuffled_edges[,"col"] <- shuffled_edges[sample(1:nrow(edges)),"col"]

# Build shuffled network
shuffled_net <- matrix(0, nrow = nrow(net), ncol = ncol(net),
                       dimnames = dimnames(net))
for(i in 1:nrow(shuffled_edges)){
  shuffled_net[shuffled_edges[i,"row"],
               shuffled_edges[i,"col"]] <- shuffled_edges[i,"Weight"]
}
rm(net,edges,shuffled_edges,i); gc()

WriteRDS(shuffled_net, paste0("Outputs/",WDIR,"/",NET,
                              "_randomized_network_Hs.rds"))


net = "stabsel"
type = "network_largest_cc"
directed = F
matrix = "Adjacency"

WDIR = "/data/public/adesous1/GeneCorrelation/"
setwd(WDIR)

.libPaths("Resources/Rlibs/R-4.0.3/")
library(igraph)

net_graph <- readRDS(paste0("Outputs/Human_Network/",net,"/",type,"_igraph.rds"))


# Get adjacency ---------------------------------------------------------------

if(matrix == "Adjacency"){
  
  if(!directed){
    # remove directionality
    net_graph <- as.undirected(net_graph, mode = "collapse", edge.attr.comb = "sum")
  }
  
  adj <- as_adjacency_matrix(net_graph, type = "both", attr = "weight",
                             sparse = FALSE)
  
  # normalize for degree / sum of weights, row-wise
  nf <- rowSums(abs(adj))
  nf[nf == 0] <- 1 # avoids dividing by 0, and keeps 0-only rows as they were
  W <- sweep(adj, 1, nf, "/")
  rm(adj,nf); gc()
  
} else if(matrix == "Laplacian"){
  
  W <- laplacian_matrix(net_graph, normalized = T)
  # laplacian dimensions are v x v
  if(isSymmetric(as.matrix(W))){
    message("A beautiful symmetric Laplacian is computed <3")
  } else{
    message("Your Laplacian is not symmetric :/")
  }
  
}
rm(net_graph); gc()

saveRDS(W, paste0("Outputs/Human_Network/",net,"/Modules/",matrix,"_weightsum_undirected_rownormalized.rds"),
        version = 2)

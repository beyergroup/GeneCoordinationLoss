.libPaths("Resources/Rlibs/R-4.0.3/")
library(devtools)
library(igraph)

directed = F
weights = "abs" # "none"

# Load BioNetSmooth package (from github)
install_github("beyergroup/BioNetSmooth")
library(BioNetSmooth)

# Load network
net_graph <- readRDS("Outputs/Human_Network/stabsel/network_largest_cc_igraph.rds")

# Remove directionality (optional)
if(!directed){
  net_graph <- as.undirected(net_graph, mode = "collapse",
                             edge.attr.comb = list(weight = function(x) x[which.max(abs(x))]))
}

# Compute weights (optional)
if(weights == "abs"){
  E(net_graph)$weight <- abs(E(net_graph)$weight)
} else if(weights == "none"){
  net_graph <- remove.edge.attribute(net_graph, "weight")
} else{
  message("nope")
}

# Convert to "edge list"
el <- as_edgelist(net_graph)

# Get batch-corrected log-normalized data for GTEx tissues 
tissue_files <- list.files("Outputs/GTEx/logFC",
                           pattern = "sampled_meanLFC.rds", full.names = T)
LFC <- sapply(tissue_files, readRDS)
colnames(LFC) <- sapply(tissue_files, function(x) strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1])
LFC <- data.frame("gene_id" = rownames(LFC), LFC)

# map LFCs to network
if(weights == "none"){
  netMap <- network_mapping(network = el, expr_mat = LFC,
                            merge.by = "gene_id",
                            global = TRUE)
} else{
  netMap <- network_mapping_weighted(network = el, expr_mat = LFC,
                                     weights = E(net_graph)$weight,
                                     merge.by = "gene_id",
                                     global = TRUE)
}

# correct network to remove topology bias effects - row-normalized adjacency
if(weights == "none"){
  adj <- as_adjacency_matrix(net_graph, type = "both",
                             sparse = FALSE)
} else{
  adj <- as_adjacency_matrix(net_graph, type = "both",
                             attr = "weight", sparse = FALSE)
}
# normalize for degree / sum of weights, row-wise
nf <- rowSums(abs(adj))
nf[nf == 0] <- 1 # avoids dividing by 0, and keeps 0-only rows as they were
W <- sweep(adj, 1, nf, "/")
rm(adj,nf); gc()

netMap$G <- W[netMap$gene_names, netMap$gene_names]
rm(W); gc()

# iterate over alphas
smoothed <- list()
for(alpha in seq(from = 0, to = 1, by = .1)){
  message("Alpha = ",alpha)
  smoothed[[as.character(alpha)]] <- network_smoothing(net = netMap$G,
                                                       mat_intensities = netMap$mat_intensities,
                                                       conditions = colnames(LFC)[-1],
                                                       iter = 50, alpha = alpha)
}
smoothed <- lapply(smoothed, function(m){
  rownames(m) <- netMap$gene_names
  return(m)})

for(alpha in names(smoothed)){
  saveRDS(smoothed[[alpha]],
          paste0("Outputs/Human_Network/stabsel/Modules/Module_Activity/Adjacency_weight",
                 weights,"_",c("undirected","directed")[directed+1],
                 "_RWRsmoothed_",alpha,"_GTEx.rds"))
}


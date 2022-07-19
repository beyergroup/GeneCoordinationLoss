# Spectral clustering on the graph matrix

args <- commandArgs(trailingOnly = TRUE)
net = args[1]
type = args[2]
weights = args[3]
matrix = args[4]

WDIR = "/data/public/adesous1/GeneCorrelation/"

setwd(WDIR)

# netfile = paste0("Outputs/Human_Network/",net,"/network_igraph.rds")
netfile = paste0("Outputs/Human_Network/",net,"/",type,"_igraph.rds")


.libPaths("Resources/Rlibs/R-4.0.3/")
library(igraph)
library(pheatmap)

# read in network
net_graph <- readRDS(netfile)

if(weights == "abs"){
  # set weights to their abs value
  net_graph <- set.edge.attribute(net_graph, "weight",
                                  value = abs(get.edge.attribute(net_graph,"weight")))
} else if(weights == "none"){
  # remove weights
  net_graph <- delete_edge_attr(net_graph, "weight")
}

# remove directionality
net_graph <- as.undirected(net_graph, mode = "collapse", edge.attr.comb = "max")

if(matrix == "Adjacency"){
  
  # get adjacency matrix
  if(weights != "none"){
    adj <- as_adjacency_matrix(net_graph, type = "both", attr = "weight", sparse = FALSE)
  } else{
    adj <- as_adjacency_matrix(net_graph, type = "both", sparse = FALSE)
  }
  
  # normalize for degree / sum of weights, row-wise
  nf <- rowSums(adj)
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

# Compute eigenvectors and eigenvalues of matrix
ee <- eigen(W, symmetric = TRUE)
rownames(ee$vectors) <- rownames(W)
saveRDS(ee, paste0("Outputs/Human_Network/",net,"/Modules/",matrix,
                   "_weight",weights,"_evectors_evalues_",type,".rds"))
min(abs(ee$values))
sum(ee$values == min(abs(ee$values)))

eigenvalues <- sort(ee$values, decreasing = FALSE)
pdf(paste0("Plots/Human_Network/",net,"/Modules/",matrix,"_weight",weights,"_evalues.pdf"),
    height = 4)
plot(eigenvalues, main = paste0(c("stabsel" = "Indirect effects included",
                                  "stabsel_pcclasso" = "Indirect effects excluded")[net],
                                ", ", matrix, " eigenvalues"),
     sub = c("abs" = "|weights|", "none" = "no weights")[weights])
dev.off()


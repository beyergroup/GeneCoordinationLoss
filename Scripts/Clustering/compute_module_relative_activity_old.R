# Compute module relative activity

args <- commandArgs(trailingOnly = TRUE)
net = "stabsel" # args[1]
type = "network_largest_cc" # args[2]
weights = "none" # args[3]
matrix = "Adjacency" # args[4]

WDIR = "/data/public/adesous1/GeneCorrelation/"
setwd(WDIR)

.libPaths("Resources/Rlibs/R-4.0.3/")
library(igraph)
source("Scripts/functions.R")

# Get module membership -------------------------------------------------------

# read in decisions
decisions <- readRDS(paste0("Outputs/Human_Network/",net,"/Modules/",matrix,
                            "_weight",weights,"_height_eigen_decision.rds"))

# read in corresponding hc object
hc <- readRDS(paste0("Outputs/Human_Network/",net,"/Modules/",matrix,"_weight",
                     weights,"_complete_link_clustering_",
                     as.character(decisions$Eigenvectors),"_evectors.rds"))

# cut at chosen height
membership <- cutree(hc, h = decisions$Height)

# remove smaller than 10 genes
modules <- names(table(membership)[table(membership) >= 10])

rm(decisions,hc); gc()


# Get network adjacency -------------------------------------------------------

net_graph <- readRDS(paste0("Outputs/Human_Network/",net,"/",type,"_igraph.rds"))

if(matrix == "Adjacency"){
  
  adj <- as_adjacency_matrix(net_graph, type = "both", attr = "weight", sparse = FALSE)
  
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


# Compute gene weights per module ---------------------------------------------

gene_weights <- list()

for(m in modules){
  
  genes <- which(membership == m)
  subW <- W[genes,genes]
  
  # matrix of network direct + indirect effects
  inverse <- solve(diag(length(genes)) - subW)
  neteffects <- subW %*% inverse
  rm(inverse); gc()
  
  # get column sum to obtain vector of gene weights
  w <- colSums(neteffects)
  
  # save in list
  gene_weights[[m]] <- w/(sum(abs(w)))
  rm(subW,inverse,neteffects); gc()
}


# Get expression data ---------------------------------------------------------

# get batch-corrected log-normalized data for GTEx tissues 

tissue_files <- list.files("GTEx_Networks/Tissue_Networks/Outputs",
                           pattern = "sampled_data.rds", full.names = T)


activity <- data.frame(row.names = modules)

for(file in tissue_files){
  
  # read in tissue-specific GTEx subset for sample identification
  gtex <- readRDS(file)
  
  # convert to gene symbol
  gtex <- DataENSGToSymbol(gtex, remove_dup = T)
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  # compute module activity
  a <- c()
  for(module in rownames(activity)){
    genes <- names(gene_weights[[module]])
    means <- rep(0,length(genes))
    names(means) <- genes
    means[intersect(genes,rownames(gtex))] <- rowMeans(gtex[intersect(genes,rownames(gtex)),])
    a <- c(a, weighted.mean(means,gene_weights[[module]]))
    
    rm(genes,means); gc()
  }
  
  activity[[tissue]] <- a
  
  rm(gtex,tissue,a); gc()
}


# save to RDS
saveRDS(activity, paste0("Outputs/Human_Network/",net,"/Modules/Module_Activity/",
                         matrix,"_weight",weights,"_relative_activity_weights_GTEx.rds"))

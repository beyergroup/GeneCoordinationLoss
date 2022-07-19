# Compute module relative activity

args <- commandArgs(trailingOnly = TRUE)
net = "stabsel" # args[1]
type = "network_largest_cc" # args[2]
weights = "sum_undirected" # args[3]
matrix = "Adjacency" # args[4]
directed = F

WDIR = "/data/public/adesous1/GeneCorrelation/"
setwd(WDIR)

.libPaths("Resources/Rlibs/R-4.0.3/")
library(igraph)
source("Scripts/functions.R")

# Get module membership -------------------------------------------------------

membership <- readRDS("Outputs/Human_Network/stabsel/Modules/old/Adjacency_weightnone_membership.rds")
modules <- names(table(membership)[table(membership) >= 10])


# Get network adjacency -------------------------------------------------------

W <- readRDS("Outputs/Human_Network/stabsel/Networks/Adjacency_undirected_weightssignedmaxabs_rownormalized_largest_cc.rds")


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

tissue_files <- list.files("Outputs/GTEx/logFC",
                           pattern = "sampled_meanLFC.rds", full.names = T)

activity <- data.frame(row.names = modules)

for(file in tissue_files){
  
  LFC <- readRDS(file)
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  # compute module activity
  a <- c()
  for(module in rownames(activity)){
    genes <- names(gene_weights[[module]])
    means <- rep(0,length(genes))
    names(means) <- genes
    means[intersect(genes,names(LFC))] <- LFC[intersect(genes,names(LFC))]
    # a <- c(a, weighted.mean(means, gene_weights[[module]]))
    a <- c(a, mean(means))
    rm(genes,means); gc()
  }
  
  activity[[tissue]] <- a
  
  rm(gtex,tissue,a); gc()
}


# save to RDS
# saveRDS(activity, paste0("Outputs/Human_Network/stabsel/Modules/Module_Activity/LocalWeightedActivity/Adjacency_undirected_weightssignedmaxabs_rownormalized_largestcc_relative_activity_weights_GTEx.rds"))
saveRDS(activity, paste0("Outputs/Human_Network/stabsel/Modules/Module_Activity/NONET_relative_activity_GTEx.rds"))

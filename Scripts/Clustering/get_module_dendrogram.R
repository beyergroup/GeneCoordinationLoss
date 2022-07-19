# Get module dendrogram from network topology

args <- commandArgs(trailingOnly = TRUE)
net = "stabsel" # args[1]
type = "network_largest_cc" # args[2]
weights = "none" # args[3]
matrix = "Adjacency" # args[4]

WDIR = "/data/public/adesous1/GeneCorrelation/"
setwd(WDIR)

.libPaths("Resources/Rlibs/R-4.0.3/")
library(dendextend)

# Get module membership -------------------------------------------------------

# read in decisions
decisions <- readRDS(paste0("Outputs/Human_Network/",net,"/Modules/",matrix,
                            "_weight",weights,"_height_eigen_decision.rds"))

# read in corresponding hc object
hc <- readRDS(paste0("Outputs/Human_Network/",net,"/Modules/",matrix,"_weight",
                     weights,"_complete_link_clustering_",
                     as.character(decisions$Eigenvectors),"_evectors.rds"))

# cut at chosen height
membership <- cutree(hc, h = decisions$Height, order_clusters_as_data = FALSE)

# remove smaller than 10 genes
modules <- names(table(membership)[table(membership) >= 10])

dend1 <- as.dendrogram(hc)
# cut
dend2 <- cut(dend1, h = decisions$Height)
# remove modules with less than 10 genes
branches <- which(unlist(lapply(dend2$lower, nleaves)) >= 10)
dend <- prune(dend2$upper, labels(dend2$upper)[-branches], reindex_dend = FALSE)

saveRDS(dend, paste0("Outputs/Human_Network/",net,"/Modules/",matrix,
                     "_weight",weights,"_pruned_dendrogram.rds"))

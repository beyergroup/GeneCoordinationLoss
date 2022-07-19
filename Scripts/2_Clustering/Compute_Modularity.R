# Compute modularity across different tree heights and eigenvector choices

args <- commandArgs(trailingOnly = TRUE)
NET = args[1]
WDIR = args[2]
PATTERN = args[3]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(igraph)
source("Scripts/functions.R")

if(PATTERN == "all"){
  files <- list.files(paste0("Outputs/",WDIR))
} else{
  # add error for when PATTERN is not among files
  files <- list.files(paste0("Outputs/",WDIR),
                      pattern = PATTERN)
  files <- files[grep(NET,files)]
}
files <- files[grep("_complete_link_clustering_",files)]

representations <- sapply(files,
                          function(x) strsplit(strsplit(x,
                                                        "_evectors_evalues")[[1]][1],
                                               paste0(NET,"_"))[[1]][2])

for(representation in unique(representations)){
  
  # get network corresponding to representation used
  net_graph <- ReadRDS(paste0("Outputs/0_Preprocessing/Adjacency/",NET,"_",
                              representation,".rds"))
  net_graph <- graph_from_adjacency_matrix(net_graph,
                  mode = c("undirected","directed")[1+GetParamFromName(representation)$Directed],
                  weighted = switch(GetParamFromName(representation)$Weights,
                                    "none" = NULL,
                                    "maxabs" = T,
                                    "signedmaxabs" = T,
                                    "sum" = T))
  
  # read in clustering results
  hc_files <- paste0("Outputs/",WDIR,"/",files[grep(representation,files)])
  hc_list <- sapply(hc_files, readRDS, simplify = F)
  names(hc_list) <- sapply(names(hc_list),
                           function(f) strsplit(f,"_")[[1]][grep("clustering",
                                                                 strsplit(f,"_")[[1]])+1])
  
  # create thresholds (100 from min to max tree heights)
  thresholds <- seq(from = 0, to = max(unlist(lapply(hc_list,
                                                     function(hc) max(hc$height)))),
                    length.out = 100)
  
  # compute modularity over thresholds and eigenvector numbers
  modularity <- lapply(hc_list, function(hc, thresholds, weights){
    clusters <- cutree(tree = hc, h = thresholds)
    if(GetParamFromName(representation)$Weights != "none"){
      modularity_scores <- apply(clusters, 2,
                                 function(m) modularity(net_graph, m,
                                                        weights = abs(E(net_graph)$weight)))
    } else {
      modularity_scores <- apply(clusters, 2,
                                 function(m) modularity(net_graph, m,
                                                        weights = NULL))
    }
    return(modularity_scores)
  }, thresholds, weights)
  modularity <- do.call(cbind,modularity)
  modularity <- modularity[,order(as.numeric(colnames(modularity)))]
  rownames(modularity) <- as.character(round(as.numeric(rownames(modularity)),2))
  
  # save whole matrix
  WriteRDS(modularity, paste0("Outputs/",WDIR,"/",NET,"_",representation,"_modularity.rds"))
  
  # find best thresholds
  cut_height = as.numeric(rownames(modularity)[which(modularity == max(modularity),
                                                     arr.ind = T)[1]])
  eigen_n = as.numeric(colnames(modularity)[which(modularity == max(modularity),
                                                  arr.ind = T)[2]])
  WriteRDS(list("Height" = cut_height, "Eigenvectors" = eigen_n),
           paste0("Outputs/",WDIR,"/",NET,"_",representation,"_height_eigen_decision.rds"))
  
  rm(net_graph, hc_files, hc_list, thresholds,
     modularity, cut_height, eigen_n); gc()
}

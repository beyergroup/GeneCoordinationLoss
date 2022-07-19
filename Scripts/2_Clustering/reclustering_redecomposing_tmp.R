# Reclustering with recomputed spectrum

NET = args[1]
WDIR = args[2]
REPRESENTATION = args[3]
THRESHOLD = args[4]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(igraph)
library(dendextend)
library(topGO)
source("Scripts/functions.R")

# Read in decision
files <- list.files(paste0("Outputs/",WDIR,"/"),
                    pattern = REPRESENTATION)
file <- files[grep("_membership",files)]
rm(files)

# Read in module calling
membership <- ReadRDS(paste0("Outputs/",WDIR,"/",NET,"_",REPRESENTATION,"_membership.rds"))

# Read in adjacency
W <- ReadRDS(paste0("Outputs/0_Preprocessing/Adjacency/",NET,"_",REPRESENTATION,".rds"))

# define eigenvector breaks
breaks <- c(5, seq(from = 10, to = 100, by = 10),
            seq(from = 150, to = 500, by = 50),
            seq(from = 600, to = 1000, by = 100),
            seq(from = 1500, to = 5000, by = 500),
            seq(from = 6000, to = 10000, by = 1000))

# problem_modules <- list("100" = c("1.1.2","5.1","8.1","1.2.3","13.1","1.2.2.1.2","1.2.4"),
#                         # "250" = c("1.1.2","5.1","8.1","1.2.3"))
#                         "250" = c("1.2.3.1.3", "1.2.2.1.2.1"))
problem_modules <- c()

while(max(table(membership[!(membership %in%
                             problem_modules)])) > THRESHOLD){
  
  # start with biggest module
  module = names(sort(table(membership[!(membership %in% problem_modules)]),
                      decreasing = T))[1]
  genes <- names(which(membership == module))
  message("Splitting module ",module," (",length(genes)," genes)")
  
  # subset to largest module
  current_W <- W[genes,genes]
  
  # Compute eigenvectors and eigenvalues of matrix
  if(grepl("undirected",file)){
    ee <- eigen(current_W, symmetric = TRUE)
  } else{
    ee <- eigen(current_W, symmetric = FALSE)
  }
  rownames(ee$vectors) <- rownames(current_W)
  
  # Hierarchical clustering
  N <- c(breaks[breaks < ncol(ee$vectors)], ncol(ee$vectors))
  # N = N[1:(length(N)-1)]
  
  hc.list <- list()
  for(n in N){
    message("Clustering on first ",n," eigenvectors")
    
    if(grepl("Adjacency",file)){
      # order based on highest absolute value
      ee_ind <- order(abs(ee$values), decreasing = T)[1:n]
    } else if(grepl("Laplacian",file)){
      # take lowest values from Laplacian (no eigenvalues below 0)
      ee_ind <- (1+(ncol(ee$vectors)-n)):ncol(ee$vectors)
    }
    d <- euc_dist(ee$vectors[,ee_ind])
    gc()
    
    if(sum(is.na(d)) > 0){
      
      # if there are NAs, use integrated dist function to re-compute those distances
      na_ind <- which(is.na(d), arr.ind = TRUE)
      dists <- apply(na_ind, 1,
                     function(x) dist(x = ee$vectors[c(x[1],x[2]),
                                                     ee_ind]))
      for(i in 1:nrow(na_ind))
        d[na_ind[i,1],na_ind[i,2]] <- dists[i]
    }
    
    dd <- as.dist(d)
    rm(d); gc()
    
    message("Distances computed")
    
    # Clustering (complete linkage)
    hc.list[[as.character(n)]] <- hclust(dd, method = "complete")
    # hc.list[[as.character(n)]] <- hclust(dd, method = "average")
    
    message("Clustering finished")
    
    rm(dd); gc()
  }
  
  # Compute modularity across hc decisions
  net_graph <- graph_from_adjacency_matrix(current_W,
                                           mode = c("undirected","directed")[1+GetParamFromName(REPRESENTATION)$Directed],
                                           weighted = switch(GetParamFromName(REPRESENTATION)$Weights,
                                                             "none" = NULL,
                                                             "maxabs" = T,
                                                             "signedmaxabs" = T,
                                                             "sum" = T))
  # handle case of no edges
  if(length(E(net_graph)) == 0){
    message("No edges left in current module: reverting to previous iteration")
    prev_module <- paste(strsplit(module,"\\.")[[1]][1:(length(strsplit(module,"\\.")[[1]])-1)],
                         collapse = ".")
    membership[sapply(membership,
                      function(x) paste(strsplit(x, "\\.")[[1]][1:length(strsplit(prev_module,"\\.")[[1]])],
                                        collapse = ".")) == prev_module] <- prev_module
    problem_modules <- c(problem_modules, prev_module)
    rm(prev_module)
    next
  }
  
  # Create tree height thresholds (100 from min to max tree heights)
  thresholds <- seq(from = 0, to = max(unlist(lapply(hc.list,
                                                     function(hc) max(hc$height)))),
                    length.out = 100)
  
  # Compute modularity
  modularity <- lapply(hc.list, function(hc, thresholds){
    clusters <- cutree(tree = hc, h = thresholds)
    if(GetParamFromName(REPRESENTATION)$Weights != "none"){
      modularity_scores <- apply(clusters, 2,
                                 function(m) modularity(net_graph, m,
                                                        weights = abs(E(net_graph)$weight)))
    } else {
      modularity_scores <- apply(clusters, 2,
                                 function(m) modularity(net_graph, m,
                                                        weights = NULL))
    }
    return(modularity_scores)
  }, thresholds)
  modularity <- do.call(cbind,modularity)
  modularity <- modularity[,order(as.numeric(colnames(modularity)))]
  rownames(modularity) <- as.character(round(as.numeric(rownames(modularity)),2))
  
  # Select optimal
  cut_height = as.numeric(rownames(modularity)[which(modularity == max(modularity),
                                                     arr.ind = T)[1,1]])
  eigen_n = colnames(modularity)[which(modularity == max(modularity),
                                       arr.ind = T)[1,2]]
  
  if(cut_height %in% c(head(thresholds,1),tail(thresholds,1))){
    message("Considering ",eigen_n," eigenvectors, cutting at h = ",cut_height)
    message("No improvement of modularity with cut: skipping module")
    problem_modules <- c(problem_modules, module)
  } else{
    message("Considering ",eigen_n," eigenvectors, cutting at h = ",cut_height)
    new_membership <- cutree(hc.list[[eigen_n]], h = cut_height)
    head(sort(table(new_membership), decreasing = T))
    
    if(length(unique(new_membership)) < 2){
      message("No new modules obtained with optimal cut: skipping module")
      problem_modules <- c(problem_modules, module)
    } else{
      # Update membership
      membership[genes] <- paste(membership[genes], new_membership[genes], sep = ".")
    }
  }
  
  message("Iteration over\n")
  
  rm(module, genes, current_W, ee, hc.list, thresholds,
     modularity, cut_height, eigen_n, new_membership)
  gc()
}


saveRDS(membership, paste0("Outputs/",WDIR,"/",NET,"_",REPRESENTATION,
                           "_membership_reclustering_redecomposing_",THRESHOLD,"_adjusted.rds"))

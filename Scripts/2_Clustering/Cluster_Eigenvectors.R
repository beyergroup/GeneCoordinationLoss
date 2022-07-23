# Hierarchical clustering of eigenvectors

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
files <- files[grep("_evectors_evalues.rds",files)]



for(file in files){
  
  # Read in eigenvalues and eigenvectors
  ee <- readRDS(paste0("Outputs/",WDIR,"/",file))
  
  N = c(5,10,20,50,100,500,1000,2000,3000,4000,5000,
        6000,7000,8000,9000,10000,15000,ncol(ee$vectors))
  
  # base dist implementation is sooooo slow, use this instead:
  euc_dist <- function(m) {
    mtm <- Matrix::tcrossprod(m)
    sq <- rowSums(m*m)
    return(sqrt(outer(sq,sq,"+") - 2*mtm))
  }
  
  for(n in N[-1]){
    
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
    
    for(met in c("complete","average")){
      hc <- hclust(dd, method = met)
      
      pdf(paste0("Plots/",WDIR,"/",gsub(".rds",
                                        paste0("_hclust_",n,"_",met,".pdf"),
                                        file)),
                 width = 20)
      plot(hc, labels = FALSE,
           main = GenerateNetworkTitle(GetParamFromName(file)),
           sub = paste0(n," eigenvectors, ",met," linkage"))
      dev.off()
      
      saveRDS(hc, paste0("Outputs/",WDIR,"/",
                         gsub(".rds",
                              paste0("_",met,"_link_clustering_",
                                     n,"_evectors.rds"),
                              file)))
      rm(hc); gc()
    }
    
    message("Clustering finished")
    
    rm(dd); gc()
  }
}

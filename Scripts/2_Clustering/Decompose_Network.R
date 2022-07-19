# Spectral clustering on the graph matrix

args <- commandArgs(trailingOnly = TRUE)
NET = args[1]
WDIR = args[2]
PATTERN = args[3]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(igraph)
source("Scripts/functions.R")

if(PATTERN == "all"){
  files <- list.files("Outputs/0_Preprocessing/Adjacency")
} else{
  # add error for when PATTERN is not among files
  files <- list.files("Outputs/0_Preprocessing/Adjacency", pattern = PATTERN)
}

for(file in files){
  
  W <- ReadRDS(paste0("Outputs/0_Preprocessing/Adjacency/",file))
  
  # Compute eigenvectors and eigenvalues of matrix
  if(grepl("undirected",file)){
    ee <- eigen(W, symmetric = TRUE)
  } else{
    ee <- eigen(W, symmetric = FALSE)
  }
  
  rownames(ee$vectors) <- rownames(W)
  rm(W); gc()
  
  # Save to file
  WriteRDS(ee, paste0("Outputs/",WDIR,"/",
                      gsub(".rds","_evectors_evalues.rds",file)))
  
  message("Min |eigenvalue| = ", min(abs(ee$values)))
  
  # Plot eigenvalues
  message("Creating eigenvalue plot")
  pdf(paste0("Plots/",WDIR,"/",gsub(".rds","_evectors_evalues.pdf",file)),
      height = 4)
  plot(ee$values,
       main = GenerateNetworkTitle(GetParamFromName(file)),
       ylab = "Ordered eigenvalues (ascending)")
  dev.off()
  
  # Plot eigengaps
  eigenvalues <- sort(abs(ee$values), decreasing = T)
  eigengap <- sapply(1:(length(eigenvalues)-1), function(i) eigenvalues[i]-eigenvalues[i+1])
  pdf(paste0("Plots/",WDIR,"/",gsub(".rds","_eigengap.pdf",file)),
      height = 4)
  plot(eigengap, main = GenerateNetworkTitle(GetParamFromName(gsub(".rds","",file))),
       ylab = "Difference to next eigenvalue",
       xlab = "Eigenvalue rank (descending order)")
  dev.off()
  
  rm(ee,eigenvalues,eigengap); gc()
}

args <- commandArgs(trailingOnly = TRUE)
NET = args[1]
WDIR = args[2]
PATTERN = args[3]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(pheatmap)
library(dendextend)
source("Scripts/functions.R")

if(PATTERN == "all"){
  files <- list.files(paste0("Outputs/",WDIR))
} else{
  # add error for when PATTERN is not among files
  files <- list.files(paste0("Outputs/",WDIR),
                      pattern = PATTERN)
  files <- files[grep(NET,files)]
}
files <- files[grep("_membership.rds",files)]

for(file in files){
  
  # Module sizes
  
  membership <- ReadRDS(paste0("Outputs/",WDIR,"/",file))
  
  sizes <- as.numeric(table(membership))
  names(sizes) <- names(table(membership))
  sizes <- sizes[sizes >= 10]
  sizes <- sort(sizes, decreasing = T)
  
  pdf(paste0("Plots/",WDIR,"/",
             gsub("membership.rds","module_sizes.pdf",file)),
      height = 5, width = 5)
  print(barplot(sizes, xlab = "Module ID", ylab = "Module size",
                main = GenerateNetworkTitle(GetParamFromName(strsplit(strsplit(file,
                                                                               paste0(NET,"_"))[[1]][2],
                                                                      "_membership")[[1]][1]))))
  dev.off()
  
  rm(sizes); gc()
  
  # Hierarchical clustering dendrogram
  
  decision <- ReadRDS(paste0("Outputs/",WDIR,"/",
                              gsub("_membership","_height_eigen_decision",file)))
  hc <- ReadRDS(paste0("Outputs/",WDIR,"/",
                       gsub("membership",
                            paste0("evectors_evalues_complete_link_clustering_",
                                   as.character(decision$Eigenvectors),
                                   "_evectors"),file)))
  
  membership <- cutree(tree = hc, h = 0.1)                                
  modules <- names(table(membership)[table(membership) >= 10])
  
  dend <- as.dendrogram(hc)
  dend_cut <- cut(dend, h = decision$Height)
  branches <- which(unlist(lapply(dend_cut$lower, nleaves)) >= 10)
  dend_pruned <- prune(dend_cut$upper, labels(dend_cut$upper)[-branches],
                       reindex_dend = FALSE)
  
  WriteRDS(dend_pruned, paste0("Outputs/",WDIR,"/",
                               gsub("membership","pruned_dendrogram",file)))
  
  pdf(paste0("Plots/",WDIR,"/",
             gsub("membership.rds","pruned_dendrogram.pdf",file)),
      width = 10, height = 5)
  print(plot(dend_pruned, leaflab = "none",
             main = GenerateNetworkTitle(GetParamFromName(strsplit(strsplit(file,
                      paste0(NET,"_"))[[1]][2],"_membership")[[1]][1]))))
  dev.off()
  
}

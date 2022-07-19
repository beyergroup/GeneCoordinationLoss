# Cluster modules based on tissue-specificity of expression

for(NET in c("stabsel","stabsel_pcclasso")){
  
  for(METHOD in c("greedy","rwalk")){
    
    m <- readRDS(list.files(paste0("Outputs/Human_Network/",NET,"/Topology/Modules"),
                            full.names = T, recursive = T, pattern = paste("mean_expression",METHOD, sep = "_")))
    if(any(is.na(m)))
      m <- m[-unique(which(is.na(m), arr.ind = T)[,1]),]
    # dd <- dist(m, method = "euclidean")
    dd <- as.dist(1-cor(t(m)))
    hc <- hclust(dd, method = "complete")
    saveRDS(hc, paste0("Outputs/Human_Network/",NET,
                       "/Topology/Modules/TissueDE/module_clustering_bytissueDE_",METHOD,".rds"))
  }
}


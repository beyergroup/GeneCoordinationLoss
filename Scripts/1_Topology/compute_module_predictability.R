# Compute predictability in each tissue + cross tissue
# (actually just subset the already computed predictability per gene)

source("Scripts/functions.R")


for(NET in c("stabsel","stabsel_pcclasso")){
  
  for(METHOD in c("greedy","rwalk")){
    
    # module membership
    membership <- readRDS(paste0("Outputs/Human_Network/",NET,
                                 "/Topology/Modules/membership_absweights_",
                                 METHOD,"_iterreclustering.rds"))
    modules <- names(table(membership)[table(membership) >= 10])
    
    # read in predictability per gene in each tissue + cross tissue
    predictability <- readRDS(paste0("Outputs/Human_Network/",NET,
                                     "/Predictability/Tissue/correlations_all.rds"))
    
    means <- matrix(nrow = length(modules), ncol = ncol(predictability),
                    dimnames = list(modules,colnames(predictability)))
    medians <- matrix(nrow = length(modules), ncol = ncol(predictability),
                      dimnames = list(modules,colnames(predictability)))
    
    for(module in modules){
      
      # restrict to genes in module
      genes <- intersect(names(which(membership == module)), rownames(predictability))
      means[module,] <- colMeans(predictability[genes,,drop = F], na.rm = T)[colnames(means)]
      medians[module,] <- apply(predictability[genes,,drop = F], 2, median, na.rm = T)[colnames(medians)]
      
    }
    
    saveRDS(means, paste0("Outputs/Human_Network/",NET,
                          "/Topology/Modules/TissuePredictability/mean_predictability_",
                          METHOD,".rds"))
    saveRDS(medians, paste0("Outputs/Human_Network/",NET,
                            "/Topology/Modules/TissuePredictability/median_predictability_",
                            METHOD,".rds"))
  }
}

# Compute within-module correlation in each tissue + cross tissue

source("Scripts/functions.R")

NET = "stabsel"
METHOD = "greedy"

# module membership
membership <- readRDS(paste0("Outputs/Human_Network/",NET,
                             "/Topology/Modules/membership_absweights_",
                             METHOD,"_iterreclustering.rds"))
modules <- names(table(membership)[table(membership) >= 10])

# log2-transformed, batch corrected, sampled GTEx data
gtex_files <- list.files("GTEx_Networks/Tissue_Networks/Outputs",
                         pattern = "sampled_data", full.names = T)
names(gtex_files) <- sapply(gtex_files, function(x)
    strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1])

means <- matrix(nrow = length(modules), ncol = length(gtex_files),
                dimnames = list(modules,names(gtex_files)))
medians <- matrix(nrow = length(modules), ncol = length(gtex_files),
                  dimnames = list(modules,names(gtex_files)))

for(tissue in names(gtex_files)){
  
  # read in tissue data
  gtex <- readRDS(gtex_files[tissue])
  gtex <- DataENSGToSymbol(data = gtex, remove_dup = T)
  
  for(module in modules){
    
    # restrict to genes in module
    genes <- names(which(membership == module))
    mat <- t(gtex[intersect(genes,rownames(gtex)),])
    
    # correlation between all genes
    cor_mat <- cor(mat, method = "pearson")
    
    # summarize
    means[module,tissue] <- mean(cor_mat[upper.tri(cor_mat)])
    medians[module,tissue] <- median(cor_mat[upper.tri(cor_mat)])
  }
}

saveRDS(means, paste0("Outputs/Human_Network/",NET,
                      "/Topology/Modules/TissueCorrelation/mean_correlation_",
                      METHOD,".rds"))
saveRDS(medians, paste0("Outputs/Human_Network/",NET,
                      "/Topology/Modules/TissueCorrelation/median_correlation_",
                      METHOD,".rds"))

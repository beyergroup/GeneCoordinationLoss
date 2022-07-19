# Compute expression in each tissue + cross tissue
# (cross tissue expression should be centered but is added
# for homogeneity with remaining analyses)

source("Scripts/functions.R")

NET = "stabsel"
METHOD = "greedy"

# module membership
membership <- readRDS(paste0("Outputs/Human_Network/",NET,
                             "/Topology/Modules/membership_absweights_",
                             METHOD,"_iterreclustering.rds"))
modules <- names(table(membership)[table(membership) >= 10])

# read in gene expression in each tissue + cross tissue
gtex_files <- list.files("GTEx_Networks/Tissue_Networks/Outputs",
                         pattern = "sampled_data.rds", full.names = T)
gtex_expr <- sapply(gtex_files, readRDS, simplify = F)
names(gtex_expr) <- sapply(names(gtex_expr),
                           function(x) strsplit(tail(strsplit(x, "/")[[1]],1),"_")[[1]][1])
gtex_expr <- lapply(gtex_expr, DataENSGToSymbol, remove_dup = T)


means <- matrix(nrow = length(modules), ncol = length(gtex_expr),
                dimnames = list(modules,names(gtex_expr)))
medians <- matrix(nrow = length(modules), ncol = length(gtex_expr),
                  dimnames = list(modules,names(gtex_expr)))
vars <- matrix(nrow = length(modules), ncol = length(gtex_expr),
               dimnames = list(modules,names(gtex_expr)))

for(module in modules){
  
  # restrict to genes in module
  genes <- names(which(membership == module))
  
  means[module,] <- unlist(lapply(gtex_expr,
    function(m) mean(m[intersect(genes,rownames(m)),])))[colnames(means)]
  medians[module,] <- unlist(lapply(gtex_expr,
    function(m) median(m[intersect(genes,rownames(m)),])))[colnames(medians)]
  vars[module,] <- unlist(lapply(gtex_expr,
    function(m) mean(apply(m[intersect(genes,rownames(m)),,drop = F], 1, var))))[colnames(vars)]
}

saveRDS(means, paste0("Outputs/Human_Network/",NET,
                      "/Topology/Modules/TissueDE/mean_expression_",
                      METHOD,".rds"))
saveRDS(medians, paste0("Outputs/Human_Network/",NET,
                        "/Topology/Modules/TissueDE/median_expression_",
                        METHOD,".rds"))
saveRDS(vars, paste0("Outputs/Human_Network/",NET,
                     "/Topology/Modules/TissueDE/mean_variance_",
                     METHOD,".rds"))

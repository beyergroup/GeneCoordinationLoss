# Dynamic tree cut

NET = "stabsel_filtered_largestCC"
WDIR = "2_Clustering"
REPRESENTATION = "Adjacency_undirected_weightssum_rownormalized"
EVECTORS = "1000"

.libPaths("Resources/Rlibs/R-4.0.3/")
library(ggplot2)
library(dynamicTreeCut)
source("Scripts/functions.R")

# Prep ------------------------------------------------------------------------

# get hierchical clustering result
hc <- ReadRDS(paste0("Outputs/",WDIR,"/",NET,"_",REPRESENTATION,
                     "_evectors_evalues_complete_link_clustering_",
                     EVECTORS,"_evectors.rds"))

# get euclidean distance based on desired evectors
ee <- ReadRDS(paste0("Outputs/",WDIR,"/",NET,"_",
                     REPRESENTATION,"_evectors_evalues.rds"))
euc_dist <- function(m) {
  mtm <- Matrix::tcrossprod(m)
  sq <- rowSums(m*m)
  return(sqrt(outer(sq,sq,"+") - 2*mtm))
}
# order based on highest absolute value
ee_ind <- order(abs(ee$values), decreasing = T)[1:as.numeric(EVECTORS)]
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

BarPlot <- function(memb, title){
  
  s <- as.numeric(table(memb))
  names(s) <- names(table(memb))
  s <- sort(s, decreasing = T)
  
  barplot(s, xlab = "Module ID", ylab = "Module size",
          main = title)
}

GetGO <- function(memb, genes, go, method = "elim"){
  
  GO.list <- list()
  
  modules <- names(sort(table(memb),decreasing = T))
  modules <- head(modules[modules != 0],10)
  
  for(module in modules){
    foreground <- genes[memb == module]
    GO.list[[module]] <- GetGOEnrich(foreground, genes,
                                     go, enrich_cutoff = 1,
                                     algorithm = method)
  }
  return(GO.list)
}

# Tree method -----------------------------------------------------------------

# tree method, deepSplit
membership_treedeep <- cutreeDynamicTree(dendro = hc,
                                         maxTreeHeight = max(hc$height),
                                         deepSplit = T, minModuleSize = 10)
# module sizes
pdf(paste0("Outputs/",WDIR,"/tmp/module_size_treedeep.pdf"),
    height = 4, width = 6)
BarPlot(membership_treedeep[membership_treedeep != 0], "tree, deep")
dev.off()
head(sort(table(membership_treedeep[membership_treedeep != 0]), decreasing = T), 20)
length(unique(membership_treedeep[membership_treedeep != 0])) # 1093 modules
length(membership_treedeep[membership_treedeep != 0])
# GO enrichment
GObp_treedeep <- GetGO(membership_treedeep, hc$labels, "BP", method = "weight01")
GOmf_treedeep <- GetGO(membership_treedeep, hc$labels, "MF", method = "weight01")
GOcc_treedeep <- GetGO(membership_treedeep, hc$labels, "CC", method = "weight01")
# save
saveRDS(membership_treedeep, paste0("Outputs/",WDIR,"/tmp/membership_treedeep.rds"))
save(file = paste0("Outputs/",WDIR,"/tmp/GO_treedeep_weight01.RData"),
     list = c("GObp_treedeep","GOmf_treedeep","GOcc_treedeep"))
rm(membership_treedeep,GObp_treedeep,GOmf_treedeep,GOcc_treedeep); gc()


# tree method, no deepSplit
membership_tree <- cutreeDynamicTree(dendro = hc,
                                     maxTreeHeight = max(hc$height),
                                     deepSplit = F, minModuleSize = 10)
# module sizes
pdf(paste0("Outputs/",WDIR,"/tmp/module_size_tree.pdf"),
    height = 4, width = 6)
BarPlot(membership_tree[membership_tree != 0], "tree, no deep")
dev.off()
head(sort(table(membership_tree[membership_tree != 0]), decreasing = T), 20)
length(unique(membership_tree[membership_tree != 0])) # 204 modules
length(membership_tree[membership_tree != 0]) # 21925 clustered genes
# GO enrichment
GObp_tree <- GetGO(membership_tree, hc$labels, "BP", method = "weight01")
GOmf_tree <- GetGO(membership_tree, hc$labels, "MF", method = "weight01")
GOcc_tree <- GetGO(membership_tree, hc$labels, "CC", method = "weight01")
# save
saveRDS(membership_tree, paste0("Outputs/",WDIR,"/tmp/membership_tree.rds"))
save(file = paste0("Outputs/",WDIR,"/tmp/GO_tree_weight01.RData"),
     list = c("GObp_tree","GOmf_tree","GOcc_tree"))
rm(membership_tree,GObp_tree,GOmf_tree,GOcc_tree); gc()

# Dynamic method --------------------------------------------------------------

res <- sapply(1:4, function(i) cutreeHybrid(dendro = hc, distM = d,
                                            deepSplit = i,
                                            minClusterSize = 10),
              simplify = FALSE)
# module sizes
pdf(paste0("Outputs/",WDIR,"/tmp/module_size_dynamic.pdf"),
    height = 4, width = 6)
lapply(res, function(l) BarPlot(l$labels[l$labels != 0], "dynamic"))
dev.off()
lapply(res, function(l) head(sort(table(l$labels[l$labels != 0]), decreasing = T), 20))
lapply(res, function(l) length(unique(l$labels[l$labels != 0]))) # 220 clusters
lapply(res, function(l) length(l$labels[l$labels != 0])) # 21947 genes clustered
GObp_dynamic <- lapply(res, function(l) GetGO(l$labels, hc$labels,
                                              "BP", method = "weight01"))
GOmf_dynamic <- lapply(res, function(l) GetGO(l$labels, hc$labels,
                                              "MF", method = "weight01"))
GOcc_dynamic <- lapply(res, function(l) GetGO(l$labels, hc$labels,
                                              "CC", method = "weight01"))
# save
saveRDS(res, paste0("Outputs/",WDIR,"/tmp/membership_dynamic.rds"))
save(file = paste0("Outputs/",WDIR,"/tmp/GO_dynamic_weight01.RData"),
     list = c("GObp_dynamic","GOmf_dynamic","GOcc_dynamic"))
rm(membership_dynamic,GObp_dynamic,GOmf_dynamic,GOcc_dynamic); gc()


save.image(paste0("Outputs/",WDIR,"/tmp/CutTreeDynamic_tmp.RData"))
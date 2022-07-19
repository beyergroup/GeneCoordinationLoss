# Apply hierarchical clustering thresholds

net = args[1]
net = "stabsel"
type = "network_largest_cc"
weights = "abs"
matrix = "Adjacency"
WDIR = "/data/public/adesous1/GeneCorrelation/"

setwd(WDIR)

.libPaths("Resources/Rlibs/R-4.0.3/")
library(igraph)

net_graph <- readRDS(paste0("Outputs/Human_Network/",net,"/",type,"_igraph.rds"))
net_graph <- set.edge.attribute(net_graph, "weight",
                                value = abs(get.edge.attribute(net_graph,"weight")))
net_graph <- as.undirected(net_graph, mode = "collapse", edge.attr.comb = "max")


hc_files <- list.files(paste0("Outputs/Human_Network/",net,"/Modules"),
                       pattern = paste0(matrix,"_weight",weights,"_complete_link_clustering_"),
                       full.names = T)

hc_list <- sapply(hc_files, readRDS, simplify = F)
names(hc_list) <- sapply(names(hc_list),
                         function(f) strsplit(f,"_")[[1]][grep("clustering",strsplit(f,"_")[[1]])+1])
thresholds <- seq(from = 0, to = max(unlist(lapply(hc_list, function(hc) max(hc$height)))), length.out = 50)

# compute modularity for given thresholds
modularity <- lapply(hc_list, function(hc, thresholds, weights){
  clusters <- cutree(tree = hc, h = thresholds)
  if(weights == "abs"){
    modularity_scores <- apply(clusters, 2, function(m) modularity(net_graph, m, weights = E(net_graph)$weight))
  } else if(weights == "none"){
    modularity_scores <- apply(clusters, 2, function(m) modularity(net_graph, m, weights = NULL))
  }
  return(modularity_scores)
}, thresholds, weights)
modularity <- do.call(cbind,modularity)
modularity <- modularity[,order(as.numeric(colnames(modularity)))]
rownames(modularity) <- as.character(round(as.numeric(rownames(modularity)),2))

library(pheatmap)
pdf(paste0("Plots/Human_Network/",net,"/Modules/",matrix,"_weight",
           weights,"_modularity_grid.pdf"), height = 8, width = 4)
pheatmap(modularity, scale = "none",
         cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 10, cellheight = 10, main = "Modularity")
dev.off()


# compute connection enrichment at given thresholds


cut_height = as.numeric(rownames(modularity)[which(modularity == max(modularity), arr.ind = T)[1]])
eigen_n = as.numeric(colnames(modularity)[which(modularity == max(modularity), arr.ind = T)[2]])

saveRDS(list("Height" = cut_height, "Eigenvectors" = eigen_n),
        paste0("Outputs/Human_Network/",net,"/Modules/",matrix,"_weight",weights,
               "_height_eigen_decision.rds"))



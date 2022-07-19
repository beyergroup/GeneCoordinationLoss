args <- commandArgs(trailingOnly = TRUE)
NET = args[1]
WDIR = args[2]
PATTERN = args[3]
REPRESENTATION = args[4]
THRESHOLD = args[5]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(igraph)
library(dendextend)
library(topGO)
library(ggpubr)
source("Scripts/functions.R")
source(paste0("Scripts/",WDIR,"/params.R"))

# Read in decision
files <- list.files(paste0("Outputs/",WDIR,"/"),
                    pattern = PATTERN)
file <- files[grep("_membership",files)]
rm(files)


# read in clustering results
hc_files <- list.files(paste0("Outputs/",WDIR),
                       pattern = PATTERN,
                       full.names = T)
hc_files <- hc_files[grep("complete_link_clustering",hc_files)]
hc_list <- sapply(hc_files, readRDS, simplify = F)
names(hc_list) <- sapply(names(hc_list),
                         function(f) strsplit(f,"_")[[1]][grep("clustering",strsplit(f,"_")[[1]])+1])

# create thresholds (100 from min to max tree heights)
thresholds <- seq(from = 0, to = max(unlist(lapply(hc_list,
                                                   function(hc) max(hc$height)))),
                  length.out = 100)

# load module calling
membership <- ReadRDS(paste0("Outputs/",WDIR,"/",file))

# load network
net_graph <- ReadRDS(paste0("Outputs/0_Preprocessing/Adjacency/",
                            gsub("_membership","",file)))
net_graph <- graph_from_adjacency_matrix(net_graph,
                                         mode = "undirected",
                                         weighted = T)

while(max(table(membership)) > THRESHOLD){
  
  # start with biggest module
  module = names(sort(table(membership), decreasing = T))[1]
  genes <- names(which(membership == module))
  
  # subset to largest module
  current_graph <- induced_subgraph(net_graph,
                                    vids = which(V(net_graph)$name %in% genes))
  
  # compute modularity over thresholds and eigenvector numbers
  modularity <- lapply(hc_list, function(hc, thresholds){
    clusters <- cutree(tree = hc, h = thresholds)
    clusters <- clusters[genes,]
    modularity_scores <- apply(clusters, 2,
                               function(m) modularity(current_graph, m,
                                                      weights = abs(E(current_graph)$weight)))
    return(modularity_scores)
  }, thresholds)
  modularity <- do.call(cbind,modularity)
  modularity <- modularity[,order(as.numeric(colnames(modularity)))]
  rownames(modularity) <- as.character(round(as.numeric(rownames(modularity)),2))
  
  # pick decision with highest modularity
  cut_height = rownames(modularity)[which(modularity == max(modularity),
                                          arr.ind = T)[1,1]]
  eigen_n = colnames(modularity)[which(modularity == max(modularity),
                                       arr.ind = T)[1,2]]
  message("Considering ",eigen_n," eigenvectors, cutting at h = ",cut_height)
  new_membership <- cutree(hc_list[[eigen_n]],
                           h = as.numeric(cut_height))[genes]
  message("New module sizes:")
  print(head(sort(table(new_membership),decreasing = T)))
  
  # update membership
  membership[names(new_membership)] <- paste(membership[names(new_membership)],
                                             new_membership, sep = ".")
  
  # clean up
  rm(module, genes, current_graph, modularity, cut_height, eigen_n, new_membership)
  gc()
}


# recompute module sizes
sizes <- as.numeric(table(membership))
names(sizes) <- names(table(membership))
sizes <- sizes[sizes >= 10]
sizes <- sort(sizes, decreasing = T)
pdf(paste0("Plots/",WDIR,"/reclustered_module_sizes_",THRESHOLD,".pdf"),
    height = 5, width = 5)
barplot(sizes, xlab = "Module ID", ylab = "Module size",
        main = paste0("Iterative re-clustering (max size = ",THRESHOLD,")"))
dev.off()

length(sizes)

# GO enrichments
GObp.list <- list()
GOmf.list <- list()
GOcc.list <- list()

for(module in names(sizes)){
  
  genes <- names(which(membership == module))
  
  GObp.list[[module]] <- GetGOEnrich(genes, names(membership),
                                     "BP", enrich_cutoff = 1,
                                      algorithm = "weight01")
  GOmf.list[[module]] <- GetGOEnrich(genes, names(membership),
                                     "MF", enrich_cutoff = 1,
                                      algorithm = "weight01")
  GOcc.list[[module]] <- GetGOEnrich(genes, names(membership),
                                     "CC", enrich_cutoff = 1,
                                     algorithm = "weight01")
}
save(GObp.list,GOmf.list,GOcc.list,
     file = paste0("Outputs/",WDIR,"/",NET,"_",REPRESENTATION,
                   "_reclustered_GO_weight01.RData"))

for(module in head(names(GObp.list),5)){
  
  plot.list <- list()
  plot.list[[1]] <- PlotGOEnrich(GObp.list[[module]], col = COLORS[module],
                                 title = paste0("GObp in module ",module))
  plot.list[[2]] <- PlotGOEnrich(GOmf.list[[module]], col = COLORS[module],
                                 title = paste0("GOcc in module ",module))
  plot.list[[3]] <- PlotGOEnrich(GOcc.list[[module]], col = COLORS[module],
                                 title = paste0("GOmf in module ",module))
  
  heights <- c(nrow(GObp.list[[module]]),
               nrow(GOmf.list[[module]]),
               nrow(GOcc.list[[module]])) + 4.5
  
  pdf(paste0("Plots/",WDIR,"/GO_weight01_module_",module,".pdf"),
      height = PDF_DIMS[module,"Height"],
      width = PDF_DIMS[module,"Width"])
  ggarrange(plotlist = plot.list, ncol = 1, nrow = length(plot.list),
            heights = heights, align = "hv")
  dev.off()
}

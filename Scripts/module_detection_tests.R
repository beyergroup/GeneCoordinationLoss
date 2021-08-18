## Comparison of module calling approaches

setwd("../")
# net = "stabsel"
net = "stabsel_pcclasso_filter01"

netfile = paste0("Outputs/Human_Network/",net,"/network_igraph.rds")

.libPaths(c("Resources/Rlibs/R-4.0.3",.libPaths()))
library(igraph, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(topGO, lib.loc = "Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")

net_graph <- readRDS(netfile)
# net_graph_undirected <- as.undirected(net_graph, mode = "collapse", edge.attr.comb = "mean")
net_graph_undirected <- as.undirected(net_graph, mode = "collapse", edge.attr.comb = "sum")
saveRDS(as_adjacency_matrix(net_graph_undirected),
        paste0("Outputs/Human_Network/",net,"/undirected_adj_mat.rds"),
        version = 2)
abs_weights <- abs(edge_attr(net_graph_undirected, "weight"))


methods <- c("greedy" = "greedy", "rwalk" = "random walk", "eigen" = "leading eigenvector")
weight_types <- c("absweights" = "|weights|", "noweights" = "no weights")


pdf(paste0("Plots/Human_Network/",net,"/Topology/community_sizes_method_compare.pdf"))
for(method in names(methods)){
  
  for(w in names(weight_types)){
    
    ww <- switch(w,
                 "absweights" = abs_weights,
                 "noweights" = NULL)
    
    modules <- switch(method,
                      "greedy" = fastgreedy.community(net_graph_undirected, weights = ww),
                      "rwalk" = walktrap.community(net_graph_undirected, weights = ww),
                      "betw" = edge.betweenness.community(net_graph_undirected, weights = ww),
                      "eigen" = leading.eigenvector.community(net_graph_undirected, weights = ww))
    
    print(barplot(sort(sizes(modules), decreasing = T),
                  names.arg = 1:length(sizes(modules)),
                  main = paste("Size of communities found with",
                               methods[method], "algorithm,",
                               weight_types[w])))
    
    saveRDS(modules, paste0("Outputs/Human_Network/",net,
                            "/Topology/Modules/modules_",w,"_",method,".rds"))
  }
}
dev.off()

rm(ww,w,method,modules); gc()



# Find communities within large communities (reclustering)

for(method in names(methods[1:2])){
  
  # read in modules
  modules <- readRDS(paste0("Outputs/Human_Network/",net,
                            "/Topology/Modules/modules_absweights_",method,".rds"))
  
  membership_log <- membership(modules)
  
  while(sum(table(membership_log) > 100) > 0){
    
    # select big modules (> 100 nodes)
    big_modules <- names(which(table(membership_log) > 100))
    
    for(module in big_modules){
      
      # reduce to current large module
      genes <- names(which(membership_log == module))
      current_graph <- induced_subgraph(net_graph_undirected, vids = genes)
      current_abs_weights <- abs(edge_attr(current_graph, "weight"))
      
      # find communities again
      current_submodules <- switch(method,
                                   "greedy" = fastgreedy.community(current_graph, weights = current_abs_weights),
                                   "rwalk" = walktrap.community(current_graph, weights = current_abs_weights),
                                   "betw" = edge.betweenness.community(current_graph, weights = current_abs_weights))
      
      # update membership
      membership_log[genes] <- paste(membership_log[genes],
                                     membership(current_submodules)[genes],
                                     sep = ".")
      
    }
    rm(big_modules,genes); gc()
  }
  
  saveRDS(membership_log, paste0("Outputs/Human_Network/",net,
                                 "/Topology/Modules/membership_absweights_",
                                 method,"_iterreclustering.rds"))
}
rm(current_abs_weights,current_graph,current_submodules,modules,membership_log); gc()



# Look at modules

pdf(paste0("Plots/Human_Network/",net,"/Topology/community_sizes_method_compare_iterreclustering.pdf"))
for(method in names(methods[1:2])){
  
  membership <- readRDS(paste0("Outputs/Human_Network/",net,
                               "/Topology/Modules/membership_absweights_",
                               method,"_iterreclustering.rds"))
  
  print(barplot(sort(table(membership), decreasing = T),
          names.arg = 1:length(table(membership)),
          main = paste("Size of communities found with",
                       methods[method], "algorithm, |weights|",
                       "\n (iterative reclustering)")))
  
  message(length(unique(membership))," modules found")
  
  # consider only > 10
  modules <- names(which(table(membership) > 10))
  message(length(modules)," of which with > 10 genes")

  # GO enrichment
  GObp.list <- list()
  GOmf.list <- list()
  for(m in modules){
    g <- names(which(membership == m))
    GObp.list[[m]] <- GetGOEnrich(g, names(V(net_graph_undirected)), "BP", enrich_cutoff = 1, algorithm = "weight")
    GOmf.list[[m]] <- GetGOEnrich(g, names(V(net_graph_undirected)), "MF", enrich_cutoff = 1, algorithm = "weight")
  }

  # save results
  save(list = c("GObp.list","GOmf.list"),
       file = paste0("Outputs/Human_Network/",net,"/Topology/Modules/GO_allcommunities_iterreclustering_",
                     method,".RData"))


  rm(membership,modules,GObp.list,GOmf.list); gc()
}
dev.off()


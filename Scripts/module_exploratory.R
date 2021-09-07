# Compare original to re-clustered modules

net = "stabsel"
# net = "stabsel_pcclasso"

.libPaths(c("Resources/Rlibs/R-4.0.3",.libPaths()))
library(igraph, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(topGO, lib.loc = "Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")

net_graph <- readRDS(paste0("Outputs/Human_Network/",net,"/network_igraph.rds"))

methods <- c("greedy" = "greedy", "rwalk" = "random walk")

for(method in names(methods)){
  
  # load original modules
  modules <- readRDS(paste0("Outputs/Human_Network/",net,"/Topology/Modules/modules_absweights_",method,".rds"))
  original_membership <- membership(modules)
  
  # 3 biggest original modules
  big_original_modules <- head(names(sort(table(original_membership), decreasing = T)), 3)
  
  # load GO enrichments after iterative reclustering
  load(paste0("Outputs/Human_Network/",net,"/Topology/Modules/GO_allcommunities_iterreclustering_",method,".RData"))
  
  # load module memberships after iterative re-clustering
  new_membership <- readRDS(paste0("Outputs/Human_Network/",net,"/Topology/Modules/membership_absweights_",
                                   method,"_iterreclustering.rds"))
  
  for(big_original_module in big_original_modules){
    
    pdf(paste0("Plots/Human_Network/",net,"/Topology/big_community_",method,"_",big_original_module,"_GObp.pdf"),
        width = 15)
    
    # enrichment in the original module
    g <- names(which(original_membership == big_original_module))
    GObp.original <- GetGOEnrich(g, names(V(net_graph)), "BP", enrich_cutoff = 1, algorithm = "weight")
    GOmf.original <- GetGOEnrich(g, names(V(net_graph)), "MF", enrich_cutoff = 1, algorithm = "weight")
    print(PlotGOEnrich(GObp.original, col = "blue",
                       title = paste0("Original module ",big_original_module," BP enrichment")))
    print(PlotGOEnrich(GOmf.original, col = "blue",
                       title = paste0("Original module ",big_original_module," MF enrichment")))
    
    color_legend <- rainbow(length(unique(new_membership[g])))
    names(color_legend) <- unique(new_membership[g])
    
    # daughter modules
    daughter_modules <- names(which(sapply(names(GObp.list),
                                           function(x) strsplit(x, "\\.")[[1]][1]) == big_original_module))
    for(daughter_module in daughter_modules){
      print(PlotGOEnrich(GObp.list[[daughter_module]], col = color_legend[daughter_module],
                         title = paste0("Re-clustered module ",daughter_module," BP enrichment")))
      print(PlotGOEnrich(GOmf.list[[daughter_module]], col = color_legend[daughter_module],
                         title = paste0("Re-clustered module ",daughter_module," MF enrichment")))
    }
    
    dev.off()
    
    
    # output daughter module members to txt file
    sink(paste0("Outputs/Human_Network/",net,"/Topology/Modules/big_community_",method,"_",
                big_original_module,"_daughtermodules.txt"), append = T)
    for(daughter_module in daughter_modules){
      print("---")
      print(daughter_module)
      print(names(new_membership[new_membership == daughter_module]))
    }
    sink()
    
    
    # visualize big cluster colored by resulting clusters
    current_graph <- induced_subgraph(net_graph, vids = which(V(net_graph)$name %in% g))
    current_graph <- set_vertex_attr(current_graph, name = "Submodule",
                                     value = new_membership[V(current_graph)$name])
    V(current_graph)$color <- color_legend[V(current_graph)$Submodule]
    E(current_graph)$width <- abs(E(current_graph)$weight)*2
    
    pdf(paste0("Plots/Human_Network/",net,"/Topology/big_community_",method,"_",
               big_original_module,"_network_viz.pdf"))
    plot(current_graph, layout = layout_with_lgl, vertex.label = NA, vertex.size = 2, edge.width = .3,
         edge.arrow.size = 0, vertex.frame.color = NA, main = paste0("Original module ",big_original_module),
         sub = "(edge direction removed)")
    for(daughter_module in daughter_modules){
      plot(induced_subgraph(current_graph, vids = which(V(current_graph)$Submodule == daughter_module)),
           layout = layout_with_lgl, vertex.size = 5, edge.curved = .1, vertex.label.dist = 1,
           edge.arrow.size = .2, vertex.frame.color = NA, vertex.label.color = "grey50",
           vertex.label.family = "Helvetica", vertex.label.cex = 1,
           main = paste0("Re-clustered module ",daughter_module))
    }
    dev.off()
  }
}



# gene group: HLA

genes <- V(net_graph)$name[grep("HLA-", V(net_graph)$name)]
# current_graph <- induced_subgraph(net_graph, vids = which(V(net_graph)$name %in% genes))
current_graph <- subgraph.edges(net_graph, eids = which(E(net_graph) %in% E(net_graph)[.inc(genes)]))
current_graph <- set_vertex_attr(current_graph, name = "Submodule",
                                 value = new_membership[V(current_graph)$name])
color_legend <- rainbow(length(unique(V(current_graph)$Submodule)))
names(color_legend) <- unique(V(current_graph)$Submodule)
V(current_graph)$color <- "black"
V(current_graph)[V(current_graph)$name %in% genes]$color <- "red"
E(current_graph)$width <- abs(E(current_graph)$weight)*2

pdf(paste0("Plots/Human_Network/",net,"/Topology/hla_network.pdf"))
plot(current_graph, layout = layout_on_sphere, vertex.size = 5,
     edge.curved = .1, vertex.label.dist = 1,
     edge.arrow.size = .2, vertex.frame.color = NA, vertex.label.color = "grey50",
     vertex.label.family = "Helvetica", vertex.label.cex = .5,
     main = "HLA genes")
V(current_graph)$color <- color_legend[V(current_graph)$Submodule]
plot(current_graph, layout = layout_on_sphere, vertex.size = 5,
     edge.curved = .1, vertex.label.dist = 1,
     edge.arrow.size = .2, vertex.frame.color = NA, vertex.label.color = "grey50",
     vertex.label.family = "Helvetica", vertex.label.cex = .5,
     main = "HLA genes", sub = "Colored by module")
dev.off()


# gene group: RBP
genes <- V(net_graph)$name[grep("RPS", V(net_graph)$name)]
genes <- c(genes, V(net_graph)$name[grep("RPL", V(net_graph)$name)])
genes_1 <- genes[grep("MRP",genes)]
genes_2 <- genes[-grep("MRP",genes)]
genes_2 <- genes_2[-grep("K", genes_2)]
genes_2 <- genes_2[-grep("TRP", genes_2)]

current_graph <- induced_subgraph(net_graph, vids = which(V(net_graph)$name %in% genes_1))
current_graph <- set_vertex_attr(current_graph, name = "Submodule",
                                 value = new_membership[V(current_graph)$name])
color_legend <- rainbow(length(unique(V(current_graph)$Submodule)))
names(color_legend) <- unique(V(current_graph)$Submodule)
V(current_graph)$color <- color_legend[V(current_graph)$Submodule]
E(current_graph)$width <- abs(E(current_graph)$weight)*2

pdf(paste0("Plots/Human_Network/",net,"/Topology/mtribosome_network_",method,"_modules.pdf"))
plot(current_graph, layout = layout_with_fr, vertex.size = 5,
     edge.curved = .1, vertex.label.dist = 1, edge.arrow.size = .2,
     vertex.frame.color = NA, vertex.label.color = "grey50",
     vertex.label.family = "Helvetica", vertex.label.cex = .8,
     main = "Mitochondrial ribosome", sub = paste0("Colored by ",methods[method]," module"))
dev.off()


current_graph <- induced_subgraph(net_graph, vids = which(V(net_graph)$name %in% genes_2))
current_graph <- set_vertex_attr(current_graph, name = "Submodule",
                                 value = new_membership[V(current_graph)$name])
color_legend <- rainbow(length(unique(V(current_graph)$Submodule)))
names(color_legend) <- unique(V(current_graph)$Submodule)
V(current_graph)$color <- color_legend[V(current_graph)$Submodule]
E(current_graph)$width <- abs(E(current_graph)$weight)*2

pdf(paste0("Plots/Human_Network/",net,"/Topology/ribosome_network_",method,"_modules.pdf"))
plot(current_graph, layout = layout_with_fr, vertex.size = 5,
     edge.curved = .1, vertex.label.dist = 1, edge.arrow.size = .2,
     vertex.frame.color = NA, vertex.label.color = "black", 
     vertex.label.family = "Helvetica", vertex.label.cex = .8,
     main = "Ribosome", sub = paste0("Colored by ",methods[method]," module"))
dev.off()

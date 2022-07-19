# Plot network subset

NET = "stabsel_filtered_largestCC"
WDIR = "2_Clustering/tmp"
# WDIR = "2_Clustering/tmp/average_linkage"
REPRESENTATION = "Adjacency_undirected_weightssum_rownormalized"
THRESHOLD = 250

.libPaths("Resources/Rlibs/R-4.0.3/")
library(igraph)
source("Scripts/functions.R")
source("Scripts/2_Clustering/params.R")

# Get membership
# membership <- ReadRDS(paste0("Outputs/2_Clustering/",NET,"_",REPRESENTATION,"_membership.rds"))
# iterative_membership <- ReadRDS(paste0("Outputs/2_Clustering/",NET,"_",REPRESENTATION,
#                                        "_membership_reclustering_redecomposing_",THRESHOLD,".rds"))
# membership <- ReadRDS(paste0("Outputs/",WDIR,"/",NET,"_",REPRESENTATION,"_membership.rds"))
# iterative_membership <- ReadRDS(paste0("Outputs/",WDIR,"/",NET,"_",REPRESENTATION,
#                                        "_membership_reclustering_redecomposing_",THRESHOLD,".rds"))
iterative_membership <- ReadRDS("Outputs/2_Clustering/stabsel_filtered_largestCC_Adjacency_undirected_weightssum_rownormalized_membership_reclustering_redecomposing_250_adjusted.rds")
membership <- sapply(iterative_membership, function(x) strsplit(x,"\\.")[[1]][1])
iterative_membership_extended <- iterative_membership
iterative_membership_extended[sapply(iterative_membership_extended,
                                     function(x) paste(strsplit(x,"\\.")[[1]][1:6],
                                                       collapse = ".")) == "1.2.3.1.1.2"] <- "1.2.3.1.1.2"
iterative_membership_extended[sapply(iterative_membership_extended,
                                     function(x) paste(strsplit(x,"\\.")[[1]][1:6],
                                                       collapse = ".")) == "1.2.2.1.2.1"] <- "1.2.2.1.2.1"
iterative_membership_extended[sapply(iterative_membership_extended,
                                     function(x) paste(strsplit(x,"\\.")[[1]][1:6],
                                                       collapse = ".")) == "1.2.4.1.1.1"] <- "1.2.4.1.1.1"
iterative_membership_extended[sapply(iterative_membership_extended,
                                     function(x) paste(strsplit(x,"\\.")[[1]][1:5],
                                                       collapse = ".")) == "1.1.2.1.1"] <- "1.1.2.1.1"
iterative_membership_extended[sapply(iterative_membership_extended,
                                     function(x) paste(strsplit(x,"\\.")[[1]][1:3],
                                                       collapse = ".")) == "5.1.1"] <- "5.1.1"
iterative_membership_extended[sapply(iterative_membership_extended,
                                     function(x) paste(strsplit(x,"\\.")[[1]][1:2],
                                                       collapse = ".")) == "8.1"] <- "8.1"
iterative_membership_extended[sapply(iterative_membership_extended,
                                     function(x) paste(strsplit(x,"\\.")[[1]][1],
                                                       collapse = ".")) == "13"] <- "13"

# Get network
net_graph <- ReadRDS(paste0("Outputs/0_Preprocessing/Adjacency/",NET,"_",REPRESENTATION,".rds"))

# Convert to graph
net_graph <- graph_from_adjacency_matrix(net_graph,
                                         mode = "undirected",
                                         weighted = T)

pdf(paste0("Plots/",WDIR,"/core_complex_modules.pdf"))
for(i in 1:length(core.complex.list)){
  
  g <- core.complex.list[[i]]
  
  # Subset to gene set
  current_graph <- induced_subgraph(net_graph, vids = which(V(net_graph)$name %in% g))
  current_graph <- set_vertex_attr(current_graph, name = "Module",
                                   value = membership[V(current_graph)$name])
  current_graph <- set_vertex_attr(current_graph, name = "Submodule",
                                   value = iterative_membership[V(current_graph)$name])
  current_graph <- set_vertex_attr(current_graph, name = "SubmoduleExtended",
                                   value = iterative_membership_extended[V(current_graph)$name])
  E(current_graph)$width <- abs(E(current_graph)$weight)*2
  
  
  # Plot
  V(current_graph)$color <- COLORS[V(current_graph)$Module]
  plot(current_graph, layout = layout_nicely,
       vertex.frame.color = "black", edge.curved = .1, vertex.label.dist = 3,
       vertex.label.color = "grey50", vertex.label.family = "Helvetica",
       main = names(core.complex.list)[i],
       sub = "One clustering round")
  if(any(!(is.na(V(current_graph)$color))))
    legend("bottomright",
           legend = names(na.omit(COLORS[unique(V(current_graph)$Module)])),
           fill = na.omit(COLORS[unique(V(current_graph)$Module)]), bty = "n")
  
  V(current_graph)$color <- COLORS[V(current_graph)$Submodule]
  plot(current_graph, layout = layout_nicely,
       vertex.frame.color = "black", edge.curved = .1, vertex.label.dist = 3,
       vertex.label.color = "grey50", vertex.label.family = "Helvetica",
       main = names(core.complex.list)[i],
       sub = "Iterative clustering")
  if(any(!(is.na(V(current_graph)$color))))
    legend("bottomright",
           legend = names(na.omit(COLORS[unique(V(current_graph)$Submodule)])),
           fill = na.omit(COLORS[unique(V(current_graph)$Submodule)]), bty = "n")
  
  V(current_graph)$color <- COLORS[V(current_graph)$SubmoduleExtended]
  plot(current_graph, layout = layout_nicely,
       vertex.frame.color = "black", edge.curved = .1, vertex.label.dist = 3,
       vertex.label.color = "grey50", vertex.label.family = "Helvetica",
       main = names(core.complex.list)[i],
       sub = "Iterative clustering extended")
  if(any(!(is.na(V(current_graph)$color))))
    legend("bottomright",
           legend = names(na.omit(COLORS[unique(V(current_graph)$SubmoduleExtended)])),
           fill = na.omit(COLORS[unique(V(current_graph)$SubmoduleExtended)]), bty = "n")
}
dev.off()

pdf(paste0("Plots/",WDIR,"/gene_family_modules.pdf"))
for(i in (1:length(gene.family.list))){
  
  g <- gene.family.list[[i]]
  
  # Subset to gene set
  if(length(intersect(V(net_graph)$name, g)) == 0)
    next
  current_graph <- induced_subgraph(net_graph, vids = which(V(net_graph)$name %in% g))
  current_graph <- set_vertex_attr(current_graph, name = "Module",
                                   value = membership[V(current_graph)$name])
  current_graph <- set_vertex_attr(current_graph, name = "Submodule",
                                   value = iterative_membership[V(current_graph)$name])
  current_graph <- set_vertex_attr(current_graph, name = "SubmoduleExtended",
                                   value = iterative_membership_extended[V(current_graph)$name])
  E(current_graph)$width <- abs(E(current_graph)$weight)*2
  
  
  # Plot
  V(current_graph)$color <- COLORS[V(current_graph)$Module]
  plot(current_graph, layout = layout_nicely,
       vertex.frame.color = "black", edge.curved = .1, vertex.label.dist = 3,
       vertex.label.color = "grey50", vertex.label.family = "Helvetica",
       main = names(gene.family.list)[i],
       sub = "One clustering round")
  if(any(!(is.na(V(current_graph)$color))))
    legend("bottomright",
           legend = names(na.omit(COLORS[unique(V(current_graph)$Module)])),
           fill = na.omit(COLORS[unique(V(current_graph)$Module)]), bty = "n")
  
  V(current_graph)$color <- COLORS[V(current_graph)$Submodule]
  plot(current_graph, layout = layout_nicely,
       vertex.frame.color = "black", edge.curved = .1, vertex.label.dist = 3,
       vertex.label.color = "grey50", vertex.label.family = "Helvetica",
       main = names(gene.family.list)[i],
       sub = "Iterative clustering")
  if(any(!(is.na(V(current_graph)$color))))
    legend("bottomright",
           legend = names(na.omit(COLORS[unique(V(current_graph)$Submodule)])),
           fill = na.omit(COLORS[unique(V(current_graph)$Submodule)]), bty = "n")
  
  V(current_graph)$color <- COLORS[V(current_graph)$SubmoduleExtended]
  plot(current_graph, layout = layout_nicely,
       vertex.frame.color = "black", edge.curved = .1, vertex.label.dist = 3,
       vertex.label.color = "grey50", vertex.label.family = "Helvetica",
       main = names(gene.family.list)[i],
       sub = "Iterative clustering extended")
  if(any(!(is.na(V(current_graph)$color))))
    legend("bottomright",
           legend = names(na.omit(COLORS[unique(V(current_graph)$SubmoduleExtended)])),
           fill = na.omit(COLORS[unique(V(current_graph)$SubmoduleExtended)]), bty = "n")
}
dev.off()


# ---- Plot individual modules ----

MODULE = "123"

g <- names(which(membership == MODULE))

# Subset to gene set
current_graph <- induced_subgraph(net_graph, vids = which(V(net_graph)$name %in% g))
current_graph <- set_vertex_attr(current_graph, name = "Module",
                                 value = membership[V(current_graph)$name])
current_graph <- set_vertex_attr(current_graph, name = "Submodule",
                                 value = iterative_membership[V(current_graph)$name])
current_graph <- set_vertex_attr(current_graph, name = "SubmoduleExtended",
                                 value = iterative_membership_extended[V(current_graph)$name])
E(current_graph)$width <- abs(E(current_graph)$weight)*20


# Plot
pdf(paste0("Plots/",WDIR,"/module_",MODULE,".pdf"))
V(current_graph)$color <- COLORS[V(current_graph)$Module]
plot(current_graph, layout = layout_with_fr,
     vertex.frame.color = "white",
     edge.curved = .1, vertex.label.dist = 2,
     vertex.label.color = "grey50", vertex.label.family = "Helvetica",
     main = paste0("Module ",MODULE))
dev.off()

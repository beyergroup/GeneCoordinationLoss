# Module identification

source("/cellnet/GeneCorrelation/Human_Network_characterization/functions.R")
.libPaths(c("/data/public/adesous1/GeneCorrelation/Resources/Rlibs/R-4.0.3/",.libPaths()))
library(igraph)
library(topGO)

setwd("../")

net_graph <- readRDS("Outputs/Human_Network/network_igraph.rds")
# remove directionality to find communities
net_graph_undirected <- as.undirected(net_graph, mode = "collapse", edge.attr.comb = "mean")

abs_weights <- abs(edge_attr(net_graph_undirected, "weight"))
modules_absweights <- cluster_louvain(graph = net_graph_undirected, weights = abs_weights)
modules_noweights <- cluster_louvain(graph = net_graph_undirected, weights = NA)
# modules_betwenness_absweights <- cluster_edge_betweenness(net_graph_undirected, weights = abs_weights)

png("Plots/Human_Network/Topology/community_sizes_louvain_abs.png", height = 300)
barplot(sort(sizes(modules_absweights), decreasing = T), names.arg = 1:length(sizes(modules_absweights)),
        main = "Size of communities found with cluster_louvain, |weights|")
dev.off()
barplot(sort(sizes(modules_noweights), decreasing = T), names.arg = 1:length(sizes(modules_noweights)),
        main = "Size of communities found with cluster_louvain, unweighted")

membership_absweights <- membership(modules_absweights)
saveRDS(membership_absweights, "Outputs/Human_Network/Topology/community_membership_louvain_abs.rds")

# remove small communities (< 3 members)
communities <- names(sort(sizes(modules_absweights), decreasing = T))
communities <- communities[communities %in% names(which(sizes(modules_absweights) > 2))]

GObp.list <- list()
GOmf.list <- list()
for(c in communities){
  g <- names(which(membership_absweights == c))
  GObp.list[[c]] <- GetGOEnrich(g, names(V(net_graph_undirected)), "BP", enrich_cutoff = 1, algorithm = "weight")
  GOmf.list[[c]] <- GetGOEnrich(g, names(V(net_graph_undirected)), "MF", enrich_cutoff = 1, algorithm = "weight")
}
save(list = c("GObp.list","GOmf.list"), file = "Outputs/Human_Network/Topology/GO_allcommunities_weight.RData")

pdf("Plots/Human_Network/Topology/GO_top10communities.pdf")
for(i in seq_along(GObp.list)[1:10]){
  # print(PlotGOEnrich(GObp.list[[i]], col = "black", title = paste("Module",names(GObp.list)[i])))
  # print(PlotGOEnrich(GOmf.list[[i]], col = "black", title = paste("Module",names(GOmf.list)[i])))
  cat(paste("\nModule",names(GObp.list)[i],"\n"))
  cat("# BIOLOGICAL PROCESSES\n")
  print(head(GObp.list[[i]][order(GObp.list[[i]]$log2Enrichment, decreasing = T),
                            c("Term","Annotated","pval","log2Enrichment")],10))
  cat("\n# MOLECULAR FUNCTIONS\n")
  print(head(GOmf.list[[i]][order(GOmf.list[[i]]$log2Enrichment, decreasing = T),
                            c("Term","Annotated","pval","log2Enrichment")]))
}
dev.off()


module_annotation <- data.frame("ModuleNumber" = names(GObp.list),
  "ModuleSize" = as.numeric(sizes(modules_absweights)[names(GObp.list)]),
  "Function" = NA)
module_annotation$Function[1] <- "Translation & RNA processing"
module_annotation$Function[2] <- "Inflammatory response"
module_annotation$Function[4] <- "Circadian signalling"
module_annotation$Function[5] <- "Response to stimulus"
module_annotation$Function[6] <- "Metabolism"
module_annotation$Function[7] <- "Neurotransmitter activity"
module_annotation$Function[8] <- "Translation & RNA processing"
module_annotation$Function[9] <- "Response to stimulus"
saveRDS(module_annotation, "objects/annotation_top10communities.rds")

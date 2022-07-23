# Single modules

source("Scripts/functions.R")
.libPaths("Resources/Rlibs/R-4.0.3/")
library(ggpubr)
library(igraph)

cluster <- "897"
colors <- c("385" = "darkgoldenrod",
            "26" = "indianred",
            "717" = "orange",
            "897" = "orange")
color <- colors[cluster]

plots <- list()
plots[[1]] <- PlotGOEnrich(GObp.list[[cluster]], col = color, title = "GO BP enrichment")
plots[[2]] <- PlotGOEnrich(GOmf.list[[cluster]], col = color, title = "GO MF enrichment")
plots[[3]] <- PlotGOEnrich(GOcc.list[[cluster]], col = color, title = "GO CC enrichment")

heights <- c(nrow(GObp.list[[cluster]]), nrow(GOmf.list[[cluster]]), nrow(GOcc.list[[cluster]]))
if(any(heights == 0)){
  plots <- plots[-(which(heights == 0))]
  heights <- heights[-which(heights == 0)]
}
heights <- heights+3.5

pdf(paste0("Plots/Human_Network/stabsel/Modules/Single_Modules/cluster_",cluster,"_GO.pdf"),
    width = 10, height = sum(heights)/3)
ggarrange(plotlist = plots, ncol = 1, heights = heights, align = "hv")
dev.off()



net_graph <- readRDS("Outputs/Human_Network/stabsel/network_igraph.rds")
membership <- readRDS("Outputs/Human_Network/stabsel/Modules/old/Adjacency_weightnone_membership.rds")
genes <- names(membership)[which(membership == cluster)]

# network module itself
current_graph <- induced_subgraph(net_graph,
                                  vids = which(V(net_graph)$name %in% genes))
E(current_graph)$width <- abs(E(current_graph)$weight)*2
V(current_graph)$color <- color
pdf(paste0("Plots/Human_Network/stabsel/Modules/Single_Modules/cluster_",cluster,"_netviz.pdf"),
    height = 5, width = 5)
PlotNet(current_graph, title = paste("Module",cluster), layout = "lgl")
dev.off()

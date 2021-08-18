setwd("../")

.libPaths(c("Resources/Rlibs/R-4.0.3",.libPaths()))
NET_FILE = "Data/Networks/Human/stabsel_network_Hs_filtered.rds" # including indirect effects
PATH = "Human_Network/stabsel/"
# NET_FILE = "Data/Networks/Human/stabsel_pcclasso_network_Hs_filtered.rds" # excluding indirect effects
# PATH = "Human_Network/stabsel_pcclasso/"

DEFAULT_COLOR = "grey"
HUB_COLOR = "darkviolet"
BOTTLENECK_COLOR = "darkmagenta"
TF_COLOR = "mediumvioletred"
POORLY_PREDICTED_COLOR = "tomato2"

source("Scripts/functions.R")

# network as a igraph object
library(igraph)
network_matrix <- as.matrix(readRDS(NET_FILE))

# make matrix square
g <- union(colnames(network_matrix),rownames(network_matrix))
adj_mat <- matrix(0, nrow = length(g), ncol = length(g), dimnames = list(g,g))
adj_mat[rownames(network_matrix),colnames(network_matrix)] <- network_matrix
rm(network_matrix,g); gc()
net_graph <- graph_from_adjacency_matrix(t(adj_mat), mode = "directed", weighted = T)
# there are no isolated vertices
saveRDS(net_graph, paste0("Outputs/",PATH,"network_igraph.rds"))

features.data <- data.frame("Gene" = names(V(graph = net_graph)))
# In-degree
features.data$InDegree <- igraph::degree(graph = net_graph, mode = "in")[features.data$Gene]
# Out-degree
features.data$OutDegree <- igraph::degree(graph = net_graph, mode = "out")[features.data$Gene]

# Betweeness (how much a given node acts as a bridge to others in the network, based on paths in the network)
features.data$Betweenness <- igraph::centr_betw(graph = net_graph)$res

# # Closeness (how connected a given node is to others in the network)
# features.data$ClosenessAll <- centr_clo(graph = net_graph, mode = "all")$res[features.data$Gene]
# features.data$ClosenessIn <- centr_clo(graph = net_graph, mode = "in")$res[features.data$Gene]
# features.data$ClosenessOut <- centr_clo(graph = net_graph, mode = "out")$res[features.data$Gene]

# Hub genes
hubs <- features.data$Gene[head(order(features.data$OutDegree, decreasing = T), 100)]
features.data$IsHub <- features.data$Gene %in% hubs

# Bottlenecks
bottlenecks <- features.data$Gene[head(order(features.data$Betweenness, decreasing = T), 100)]
features.data$IsBottleneck <- features.data$Gene %in% bottlenecks

# TFs
TFtable <- read.delim("Resources/TFtable.txt")
TFs <- as.character(TFtable$Gene_symbol)
features.data$IsTF <- features.data$Gene %in% TFs

table(features.data$IsHub, features.data$IsTF)
table(features.data$IsBottleneck, features.data$IsTF)
table(features.data$IsHub, features.data$IsBottleneck)

saveRDS(hubs, paste0("Outputs/",PATH,"Topology/100_hubs.rds"))
saveRDS(bottlenecks, paste0("Outputs/",PATH,"/Topology/100_bottlenecks.rds"))

poorly_predicted <- readRDS(paste0("Outputs/",PATH,"Predictability/Tissue/poorly_predicted_crosstissue.rds"))
features.data$IsPoorlyPredicted <- features.data$Gene %in% poorly_predicted

library(ggplot2)
library(ggpubr)

plot.list <- list()
# All genes
plot.list[[1]] <- ggplot(features.data) +
  geom_bar(aes(x = InDegree), fill = DEFAULT_COLOR, color = DEFAULT_COLOR) +
  ggtitle("All nodes") + xlab("In-degree") +
  theme_classic() + theme(text = element_text(size = 20))
plot.list[[2]] <- ggplot(features.data) +
  geom_bar(aes(x = OutDegree), fill = DEFAULT_COLOR, color = DEFAULT_COLOR) +
  ggtitle("All nodes")+ xlab("Out-degree") +
  theme_classic() + theme(text = element_text(size = 20))
plot.list[[3]] <- ggplot(features.data) +
  geom_density(aes(x = Betweenness), fill = DEFAULT_COLOR, color = DEFAULT_COLOR) +
  xlab("Betweeness") + ggtitle("All nodes") +
  theme_classic() + theme(text = element_text(size = 20),
    axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
# Hubs
plot.list[[4]] <- ggplot(subset(features.data, IsHub)) +
  geom_bar(aes(x = InDegree), fill = HUB_COLOR, color = HUB_COLOR) +
  xlim(c(min(features.data$InDegree),max(features.data$InDegree))) +
  ggtitle("Hub nodes") + xlab("In-degree") +
  theme_classic() + theme(text = element_text(size = 20))
plot.list[[5]] <- ggplot(subset(features.data, IsHub)) +
  geom_bar(aes(x = OutDegree), fill = HUB_COLOR, color = HUB_COLOR) +
  xlim(c(min(features.data$OutDegree),max(features.data$OutDegree))) +
  xlab("Out-degree") + ggtitle("Hub nodes") +
  theme_classic() + theme(text = element_text(size = 20))
plot.list[[6]] <- ggplot(subset(features.data, IsHub)) +
  geom_density(aes(x = Betweenness), fill = HUB_COLOR, color = HUB_COLOR) +
  xlim(c(min(features.data$Betweenness),max(features.data$Betweenness))) +
  xlab("Betweenness") + ggtitle("Hub nodes") +
  theme_classic() + theme(text = element_text(size = 20),
    axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
# Bottlenecks
plot.list[[7]] <- ggplot(subset(features.data, IsBottleneck)) +
  geom_bar(aes(x = InDegree), fill = BOTTLENECK_COLOR, color = BOTTLENECK_COLOR) +
  xlim(c(min(features.data$InDegree),max(features.data$InDegree))) +
  ggtitle("Bottlenecks") + xlab("In-degree") +
  theme_classic() + theme(text = element_text(size = 20))
plot.list[[8]] <- ggplot(subset(features.data, IsBottleneck)) +
  geom_bar(aes(x = OutDegree), fill = BOTTLENECK_COLOR, color = BOTTLENECK_COLOR) +
  xlim(c(min(features.data$OutDegree),max(features.data$OutDegree))) +
  xlab("Out-degree") + ggtitle("Bottlenecks") +
  theme_classic() + theme(text = element_text(size = 20))
plot.list[[9]] <- ggplot(subset(features.data, IsBottleneck)) +
  geom_density(aes(x = Betweenness), fill = BOTTLENECK_COLOR, color = BOTTLENECK_COLOR) +
  xlim(c(min(features.data$Betweenness),max(features.data$Betweenness))) +
  xlab("Betweenness") + ggtitle("Bottlenecks") +
  theme_classic() + theme(text = element_text(size = 20),
    axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
# TFs
plot.list[[10]] <- ggplot(subset(features.data, IsTF)) +
  geom_bar(aes(x = InDegree), fill = TF_COLOR, color = TF_COLOR) +
  xlim(c(min(features.data$InDegree),max(features.data$InDegree))) +
  ggtitle("TFs") + xlab("In-degree") +
  theme_classic() + theme(text = element_text(size = 20))
plot.list[[11]] <- ggplot(subset(features.data, IsTF)) +
  geom_bar(aes(x = OutDegree), fill = TF_COLOR, color = TF_COLOR) +
  xlim(c(min(features.data$OutDegree),max(features.data$OutDegree))) +
  xlab("Out-degree") + ggtitle("TFs") +
  theme_classic() + theme(text = element_text(size = 20))
plot.list[[12]] <- ggplot(subset(features.data, IsTF)) +
  geom_density(aes(x = Betweenness), fill = TF_COLOR, color = TF_COLOR) +
  xlim(c(min(features.data$Betweenness),max(features.data$Betweenness))) +
  xlab("Betweenness") + ggtitle("TFs") +
  theme_classic() + theme(text = element_text(size = 20),
    axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
# poorly predicted
plot.list[[13]] <- ggplot(subset(features.data, IsPoorlyPredicted)) +
  geom_bar(aes(x = InDegree), fill = POORLY_PREDICTED_COLOR, color = POORLY_PREDICTED_COLOR) +
  xlim(c(min(features.data$InDegree),max(features.data$InDegree))) +
  ggtitle("Poorly predicted") + xlab("In-degree") +
  theme_classic() + theme(text = element_text(size = 20))
plot.list[[14]] <- ggplot(subset(features.data, IsPoorlyPredicted)) +
  geom_bar(aes(x = OutDegree), fill = POORLY_PREDICTED_COLOR, color = POORLY_PREDICTED_COLOR) +
  xlim(c(min(features.data$OutDegree),max(features.data$OutDegree))) +
  xlab("Out-degree") + ggtitle("Poorly predicted") +
  theme_classic() + theme(text = element_text(size = 20))
plot.list[[15]] <- ggplot(subset(features.data, IsPoorlyPredicted)) +
  geom_density(aes(x = Betweenness), fill = POORLY_PREDICTED_COLOR, color = POORLY_PREDICTED_COLOR) +
  xlim(c(min(features.data$Betweenness),max(features.data$Betweenness))) +
  xlab("Betweenness") + ggtitle("Poorly predicted") +
  theme_classic() + theme(text = element_text(size = 20),
    axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())

pdf(paste0("Plots/",PATH,"Topology/centrality_all_hubs_bottlenecks_TFs_poorly_predicted.pdf"),
    # height = 10, width = 12)
    height = 12, width = 12)
# ggarrange(plotlist = plot.list, align = "hv", ncol = 3, nrow = 4)
ggarrange(plotlist = plot.list, align = "hv", ncol = 3, nrow = 5)
dev.off()


# GO enrichments in hubs and bottlenecks

hubs_BP <- GetGOEnrich(hubs, names(V(net_graph)), "BP")
hubs_MF <- GetGOEnrich(hubs, names(V(net_graph)), "MF")
bottlenecks_BP <- GetGOEnrich(bottlenecks, names(V(net_graph)), "BP")
bottlenecks_MF <- GetGOEnrich(bottlenecks, names(V(net_graph)), "MF")
save(list = c("hubs_BP","hubs_MF","bottlenecks_BP","bottlenecks_MF"),
     file = paste0("Outputs/",PATH,"Topology/GO_hubs_bottlenecks.RData"))

go.list <- list()
go.list[[1]] <- PlotGOEnrich(hubs_BP, HUB_COLOR, "GObp enriched in hubs")
go.list[[2]] <- PlotGOEnrich(hubs_MF, HUB_COLOR, "GOmf enriched in hubs")
go.list[[3]] <- PlotGOEnrich(bottlenecks_BP, BOTTLENECK_COLOR, "GObp enriched in bottlenecks")
go.list[[4]] <- PlotGOEnrich(bottlenecks_MF, BOTTLENECK_COLOR, "GOmf enriched in bottlenecks")

if(grepl("pcclasso",PATH)){
  
  if(grepl("filter01",PATH)){
    
    pdf(paste0("Plots/",PATH,"Topology/GO_hubs.pdf"), width = 13, height = 12)
    print(ggarrange(plotlist = go.list[1:2], ncol = 1, heights = c(4,1), align = "hv"))
    dev.off()
    
    pdf(paste0("Plots/",PATH,"Topology/GO_bottlenecks.pdf"), width = 13, height = 7)
    print(ggarrange(plotlist = go.list[3:4], ncol = 1, heights = c(2.5,1), align = "hv"))
    dev.off()
    
  } else{
    
    pdf(paste0("Plots/",PATH,"Topology/GO_hubs.pdf"), width = 13, height = 12)
    ggarrange(plotlist = go.list[1:2], ncol = 1, heights = c(4,1), align = "hv")
    dev.off()
    
    pdf(paste0("Plots/",PATH,"Topology/GO_bottlenecks.pdf"), width = 13, height = 7)
    ggarrange(plotlist = go.list[3:4], ncol = 1, heights = c(2.5,1), align = "hv")
    dev.off()
    
  }
  
} else{
  
  pdf(paste0("Plots/",PATH,"Topology/GO_hubs.pdf"), width = 13, height = 8)
  ggarrange(plotlist = go.list[1:2], ncol = 1, heights = c(7,2), align = "hv")
  dev.off()
  
  pdf(paste0("Plots/",PATH,"Topology/GO_bottlenecks.pdf"), width = 13, height = 10)
  ggarrange(plotlist = go.list[3:4], ncol = 1, heights = c(3,1), align = "hv")
  dev.off()
  
}

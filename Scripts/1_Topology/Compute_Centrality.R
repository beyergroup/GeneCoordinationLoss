# Compute centrality measures

args <- commandArgs(trailingOnly = T)
NET = args[1]
WDIR = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")
library(igraph)


# Read in network and convert to igraph ---------------------------------------

net <- ReadRDS(paste0("Outputs/0_Preprocessing/",NET,"_network_Hs.rds"))
net <- MatrixToSquare(net)
g <- AdjToIgraph(net)


# Compute centrality measures -------------------------------------------------

features.data <- data.frame("Gene" = names(V(graph = g)))

# In-degree
features.data[["In-Degree"]] <- igraph::degree(graph = g, mode = "in")[features.data$Gene]
# Out-degree
features.data[["Out-Degree"]] <- igraph::degree(graph = g, mode = "out")[features.data$Gene]
# Betweenness (how much a given node acts as a bridge to others in the network,
# based on paths in the network)
features.data[["Betweenness"]] <- igraph::centr_betw(graph = g)$res


# Identify relevant genes -----------------------------------------------------

# Hub genes
hubs <- features.data$Gene[head(order(features.data[["Out-Degree"]],
                                      decreasing = T), 100)]
features.data[["Hub genes"]] <- features.data$Gene %in% hubs

# Bottlenecks
bottlenecks <- features.data$Gene[head(order(features.data[["Betweenness"]],
                                             decreasing = T), 100)]
features.data$Bottlenecks <- features.data$Gene %in% bottlenecks

# TFs
TFtable <- read.delim("Resources/TFtable.txt")
TFs <- as.character(TFtable$Gene_symbol)
features.data$TFs <- features.data$Gene %in% TFs


# Compare subsets -------------------------------------------------------------

message(sum(features.data[["Hub genes"]] & features.data$TFs),"/",
        sum(features.data[["Hub genes"]])," TFs in Hubs compared to ",
        sum(features.data$TFs),"/",nrow(features.data)," in background")

message(sum(features.data$Bottlenecks & features.data$TFs),"/",
        sum(features.data$Bottleneck)," TFs in Bottlenecks compared to ",
        sum(features.data$TFs),"/",nrow(features.data)," in background")

message(sum(features.data[["Hub genes"]] & features.data$Bottlenecks),
        " Hubs and Bottlenecks in common")


# Save outputs ----------------------------------------------------------------
WriteRDS(hubs, paste0("Outputs/",WDIR,"/100_hubs_",NET,".rds"))
WriteRDS(bottlenecks, paste0("Outputs/",WDIR,"/100_bottlenecks_",NET,".rds"))
WriteRDS(features.data, paste0("Outputs/",WDIR,"/centrality_measures_",NET,".rds"))


.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")
library(igraph)
source(paste0("Scripts/1_Topology/params.R"))
library(ggplot2)
library(ggpubr)

MODE = "cis"


# Read in network and convert to igraph ---------------------------------------

net <- ReadRDS(paste0("Outputs/5_Predictability/",MODE,"_net_chrarm.rds"))
net <- net[rowSums(net != 0) != 0, colSums(net != 0) != 0]
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
features.data[["Hub nodes"]] <- features.data$Gene %in% hubs

# Bottlenecks
bottlenecks <- features.data$Gene[head(order(features.data[["Betweenness"]],
                                             decreasing = T), 100)]
features.data$Bottlenecks <- features.data$Gene %in% bottlenecks

# TFs
TFtable <- read.delim("Resources/TFtable.txt")
TFs <- as.character(TFtable$Gene_symbol)
features.data$TFs <- features.data$Gene %in% TFs


# Compare subsets -------------------------------------------------------------

message(sum(features.data[["Hub nodes"]] & features.data$TFs),"/",
        sum(features.data[["Hub nodes"]])," TFs in Hubs compared to ",
        sum(features.data$TFs),"/",nrow(features.data)," in background")

message(sum(features.data$Bottlenecks & features.data$TFs),"/",
        sum(features.data$Bottleneck)," TFs in Bottlenecks compared to ",
        sum(features.data$TFs),"/",nrow(features.data)," in background")

message(sum(features.data[["Hub nodes"]] & features.data$Bottlenecks),
        " Hubs and Bottlenecks in common")

# Create list of plot arguments -----------------------------------------------

centrality_measures <- names(which(lapply(features.data,class) == "numeric"))
node_types <- c("All nodes",
                names(which(lapply(features.data,class) == "logical")))
args.mat <- expand.grid(x = centrality_measures,
                        Title = node_types,
                        stringsAsFactors = F)
args.mat$Color <- COLOR_SCHEME[args.mat$Title]
args.mat$Subset <- args.mat$Title != "All nodes"


# Apply plotting functions ----------------------------------------------------

plot.list <- apply(args.mat, 1,
                   function(args){
                     if(args["x"] == "Betweenness"){
                       DensityPlot(features.data, args)
                     } else{
                       BarPlot(features.data, args)
                     }})

pdf(paste0("Plots/5_Predictability/",MODE,"_centrality.pdf"),
    height = DIMS["All","Height"], width = DIMS["All","Width"])
ggarrange(plotlist = plot.list,
          ncol = length(centrality_measures),
          nrow = length(node_types))
dev.off()

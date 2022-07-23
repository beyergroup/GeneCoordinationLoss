args <- commandArgs(trailingOnly = T)
NET = args[1]
WDIR = args[2]
POORLY_PREDICTED = as.logical(args[3])

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")
source(paste0("Scripts/",WDIR,"/params.R"))
library(ggplot2)
library(ggpubr)

# Read in centrality measures -------------------------------------------------
features.data <- ReadRDS(paste0("Outputs/",WDIR,"/centrality_measures_",
                                NET,".rds"))

# Create list of plot arguments -----------------------------------------------
centrality_measures <- names(which(lapply(features.data,class) == "numeric"))
node_types <- c("All genes",
                names(which(lapply(features.data,class) == "logical")))
args.mat <- expand.grid(x = centrality_measures,
                        Title = node_types,
                        stringsAsFactors = F)
args.mat$Color <- COLOR_SCHEME[args.mat$Title]
args.mat$Subset <- args.mat$Title != "All genes"

# Apply plotting functions ----------------------------------------------------
plot.list <- apply(args.mat, 1,
                   function(args){
                     if(args["x"] == "Betweenness"){
                       DensityPlot(features.data, args)
                     } else{
                       BarPlot(features.data, args)
                     }})

# Save to pdf -----------------------------------------------------------------
pdf(paste0("Plots/",WDIR,"/centrality_",NET,
           c(".pdf","_poorlypredicted.pdf")[POORLY_PREDICTED+1]),
    height = DIMS[POORLY_PREDICTED+1,"Height"],
    width = DIMS[POORLY_PREDICTED+1,"Width"])
ggarrange(plotlist = plot.list,
          ncol = length(centrality_measures),
          nrow = length(node_types))
dev.off()

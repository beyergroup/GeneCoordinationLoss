# Run GO enrichment analysis

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

# Create list of arguments ----------------------------------------------------
gene_groups <- names(which(lapply(features.data,class) == "logical"))
gene_groups <- gene_groups[-grep("TF",gene_groups)]
if((!POORLY_PREDICTED) & ("Poorly predicted" %in% gene_groups)){
  gene_groups <- gene_groups[-grep("Poorly predicted",gene_groups)]
}
args.mat <- expand.grid(Group = gene_groups,
                        GO_type = c("BP","MF"),
                        stringsAsFactors = F)
args.mat$Color <- COLOR_SCHEME[args.mat$Group]


# GO enrichments --------------------------------------------------------------

GO.list <- apply(args.mat, 1,
                 function(args){
                   foreground <- features.data[features.data[[args[["Group"]]]],
                                               "Gene"]
                   GetGOEnrich(foreground, background = features.data$Gene,
                               go = args[["GO_type"]])})
names(GO.list) <- paste0(args.mat$Group," (",args.mat$GO_type,")")

WriteRDS(GO.list,
         file = paste0("Outputs/",WDIR,"/centrality_GO_",NET,
                       c(".rds","_poorlypredicted.rds")[POORLY_PREDICTED+1]))


# Plot ------------------------------------------------------------------------

plot.list <- sapply(1:length(GO.list),
                    function(i) PlotGOEnrich(GO.list[[i]],
                                             args.mat[i,"Color"],
                                             paste0(args.mat[i,"Group"]," (",
                                                    args.mat[i,"GO_type"],")")),
                    simplify = F)

for(g in unique(args.mat$Group)){
  index <- which(args.mat$Group == g)
  heights <- unlist(lapply(GO.list[index], nrow)) + 3.5
  pdf(paste0("Plots/",WDIR,"/centrality_GO_",NET,"_",gsub(" ","_",g),
             c(".pdf","_poorlypredicted.pdf")[POORLY_PREDICTED+1]),
      height = subset(GO_PDF_DIMS, (Net == NET) & (Group == g))$Height,
      width = subset(GO_PDF_DIMS, (Net == NET) & (Group == g))$Width)
  print(ggarrange(plotlist = plot.list[index], ncol = 1, nrow = 2,
                  heights = heights, align = "hv"))
  dev.off()
}

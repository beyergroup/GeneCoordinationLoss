# Compare in-going edges to out-going edges

args <- commandArgs(trailingOnly = T)
NET = args[1]
WDIR = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(ggplot2)
source("Scripts/functions.R")

# Load network ----------------------------------------------------------------
net <- as.matrix(ReadRDS(paste0("Outputs/",WDIR,"/",NET,"_network_Hs.rds")))

# Collect values for edge pairs (node i -> node j and opposite direction) -----
ind <- which(upper.tri(net), arr.ind = T)
e <- apply(ind, 1,
           function(x) (net[x[1],x[2]] == 0) & (net[x[2],x[1]] == 0))

# Write percentage of pairs that are non-zero in both directions ------------------
f <- file(paste0("Tables/",WDIR,"/percentnonzeropairs.txt"))
writeLines(as.character(sum(!e)/length(e)*100), f)
close(f)

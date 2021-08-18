# Investigate role of nuclear lamins as bottlenecks in the network

.libPaths(c("Resources/Rlibs/R-4.0.3",.libPaths()))
library(biomaRt)
source("Scripts/functions.R")

bottlenecks <- readRDS("/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel_pcclasso/Topology/100_bottlenecks.rds")

lamins <- GetGOGenes("GO:0005652")
lamins <- unique(lamins$hgnc_symbol)
# lamin genes are not bottlenecks (probably only lamin interactors)
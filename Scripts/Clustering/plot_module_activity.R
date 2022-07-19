# plot activity data.frames

args <- commandArgs(trailingOnly = TRUE)
net = "stabsel" # args[1]
# type = "network_largest_cc" # args[2]
weights = "abs" # "ident" # "sum_undirected" # args[3]
matrix = "Adjacency" # "directed_full" # args[4]
# type = "relative_activity" # "activity_weights"

WDIR = "/data/public/adesous1/GeneCorrelation/"
setwd(WDIR)

.libPaths("Resources/Rlibs/R-4.0.3/")
library(pheatmap)
library(RColorBrewer)

dend <- readRDS(paste0("Outputs/Human_Network/",net,"/Modules/Adjacency_weightnone_pruned_dendrogram.rds"))

for(alpha in seq(from = 0, to = 1, by = .1)){
  
  for(directed in c(TRUE,FALSE)){
    
    type = paste0(c("undirected","directed")[directed+1],
                  "_RWRsmoothed_",alpha,"_relative_activity")
    
    file <- paste0("Outputs/Human_Network/",net,"/Modules/Module_Activity/",
                  matrix,"_weight",weights,"_",type,"_GTEx.rds")
    
    activity <- readRDS(file)
    activity <- activity[,-grep("CrossTissue",colnames(activity))]
    plot.data <- as.matrix(activity)
    
    pheatmap(t(plot.data), scale = "none", border_color = NA,
             color = colorRampPalette(c("blue","white","red"))(n = 501),
             breaks = seq(from = -max(abs(plot.data)), to = max(abs(plot.data)), length.out = 501),
             cellwidth = 10, cellheight = 10, show_colnames = T,
             cluster_cols = as.hclust(dend),
             main = c("relative_activity" = "Relative module activity in GTEx tissues",
                      "activity" = "Module activity in GTEx tissues")[type],
             filename = paste0("Plots/Human_Network/stabsel/Modules/Module_Activity/",
                               matrix,"_weight",weights,"_",type,"_GTEx.pdf"))
  }
  
}


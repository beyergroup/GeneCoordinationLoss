# Plot modularity heatmaps

args <- commandArgs(trailingOnly = TRUE)
NET = args[1]
WDIR = args[2]
PATTERN = args[3]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(pheatmap)
source("Scripts/functions.R")

if(PATTERN == "all"){
  files <- list.files(paste0("Outputs/",WDIR))
} else{
  # add error for when PATTERN is not among files
  files <- list.files(paste0("Outputs/",WDIR),
                      pattern = PATTERN)
  files <- files[grep(NET,files)]
}
files <- files[grep("_modularity",files)]


for(file in files){
  
  modularity <- ReadRDS(paste0("Outputs/",WDIR,"/",file))
  
  pdf(paste0("Plots/",WDIR,"/",
             gsub("modularity.rds","modularity_grid.pdf",file)),
      height = 12, width = 4)
  print(pheatmap(modularity, scale = "none",
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 cellwidth = 10, cellheight = 8, main = "Modularity"))
  dev.off()
  
  pdf(paste0("Plots/",WDIR,"/",
             gsub("modularity.rds","modularity_trend.pdf",file)),
      height = 4, width = 14)
  plot(as.numeric(modularity),
       xlab = "Eigenvector thresholds",
       ylab = "Modularity",
       xaxt = "none", pch = 20, cex = .6)
  axis(1, cumsum(rep(nrow(modularity),
                     ncol(modularity)))-50,
       labels = colnames(modularity),
       tick = FALSE)
  abline(v = c(0,cumsum(rep(nrow(modularity),
                            ncol(modularity)))),
         lty = "dashed")
  dev.off()
  
}

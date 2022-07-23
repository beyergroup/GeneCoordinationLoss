
library(pheatmap, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(RColorBrewer, lib.loc = "Resources/Rlibs/R-4.0.3/")

method_labels <- c("greedy" = "Greedy", "rwalk" = "Random Walk")

for(NET in c("stabsel","stabsel_pcclasso")){
  
  for(METHOD in c("greedy","rwalk")){
    
    files <- list.files(paste0("Outputs/Human_Network/",NET,"/Topology/Modules/TissueDE"),
                        pattern = paste0("absexpression_",METHOD,".rds"), full.names = T)
    
    abs_expr <- sapply(files,readRDS)
    colnames(abs_expr) <- sapply(colnames(abs_expr), function(x) strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1])
    if(sum(is.na(abs_expr)) > 0)
      abs_expr <- abs_expr[-unique(which(is.na(abs_expr), arr.ind = T)[,1]),]
    
    pheatmap(t(abs_expr), scale = "none",
             color = colorRampPalette(c("white","red"))(n = 501),
             cellwidth = 1.5, cellheight = 10, show_colnames = F,
             main = paste0("Absolute expression (",
                           method_labels[METHOD]," modules)"),
             filename = paste0("Plots/Human_Network/",NET,"/Topology/Modules/absexpression_",
                               METHOD,".pdf"))
  }
}

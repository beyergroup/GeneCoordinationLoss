SLOPE_DIR = "Outputs/5_Predictability/Age"
TISSUES = "all"

library(ggplot2)
library(ggpubr)
source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")


# get slopes with 70s
slope_files <- list.files(SLOPE_DIR,
                          pattern = paste0("ageslope_well_predicted.rds"),
                          full.names = T)
slopes <- sapply(slope_files, ReadRDS, simplify = F)
names(slopes) <- sapply(names(slopes),
                        function(x) strsplit(tail(strsplit(x, "/")[[1]],1),
                                             "_")[[1]][1])
if(TISSUES == "main"){
  slope_files <- list.files(SLOPE_DIR,
                            pattern = paste0("ageslope_well_predicted_no70.rds"),
                            full.names = T)
  slopes_no70 <- sapply(slope_files, ReadRDS, simplify = F)
  names(slopes_no70) <- sapply(names(slopes_no70),
                               function(x) strsplit(tail(strsplit(x, "/")[[1]],1),
                                                    "_")[[1]][1])
}

rm(slope_files); gc()


for(tissue in names(slopes)){
  
  data <- slopes[[tissue]]
  
  if(TISSUES == "main"){
    
    data_no70 <- slopes_no70[[tissue]]
    
    # reduce to common genes
    genes <- intersect(rownames(data),rownames(data_no70))
    data <- data[genes,]
    data_no70 <- data_no70[genes,]
    
    # reduce to genes with same sign on both analyses
    genes <- names(which(sign(data[,"Slope"]) == sign(data_no70[,"Slope"])))
    data <- data[genes,]
    data_no70 <- data_no70[genes,]
    
    # get rank of each gene on both analyses according to pvalue
    ranks <- data.frame("Gene" = genes, "Full" = rank(data[,"pval"]),
                        "No70" = rank(data_no70[,"pval"]))
    
    # get top 100
    ranks <- ranks[order(rowMeans(ranks[,c("Full","No70")])),]
    
  } else{
    ranks <- data.frame("Gene" = rownames(data), "pval" = rank(data[,"pval"]))
    ranks <- ranks[order(ranks[,"pval"]),]
  }
  
  hits <- ranks[1:100,"Gene"]
  
  hits_up <- names(which(data[hits, "Slope"] > 0))
  hits_dw <- names(which(data[hits, "Slope"] < 0))
  
  if(TISSUES == "all"){
    WriteRDS(list("Up" = hits_up, "Down" = hits_dw),
             paste0(SLOPE_DIR,"/",tissue,"_predictability_hits.rds"))
  } else{
    WriteRDS(list("Up" = hits_up, "Down" = hits_dw),
             paste0(SLOPE_DIR,"/",tissue,"_predictability_hits_robust.rds"))
  }
}

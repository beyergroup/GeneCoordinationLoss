# Compute correlations between network predictions and original data

args <- commandArgs(trailingOnly = T)
EXPR_DIR = "Outputs/5_Predictability" # args[1]
PRED_DIR = "Outputs/5_Predictability" # args[2]
COR = "spearman"

source("Scripts/functions.R")

files <- list.files(PRED_DIR, pattern = "net_predictions.rds")
files <- files[grep("young|old",files)]
# files <- files[-grep("young|old",files)]
data_files <- list.files(EXPR_DIR, pattern = "centered_data.rds")
# data_files <- data_files[grep("young|old",data_files)]

for(mode in c("cis","trans")){
  
  f <- files[grep(mode,files)]
  
  for(i in 1:length(f)){
    
    predicted <- as.matrix(readRDS(paste0(PRED_DIR,"/",f[i])))
    centered_data <- as.matrix(readRDS(paste0(EXPR_DIR,"/",data_files[i])))
    centered_data <- DataENSGToSymbol(centered_data)
    centered_data <- centered_data[rownames(predicted),]
    
    if(sum(is.na(predicted)) > 0)
      break
    
    cor <- switch(COR,
                  "pearson" = sapply(rownames(predicted),
                                     function(gene_name) cor.test(x = centered_data[gene_name,],
                                                                  y = predicted[gene_name,],
                                                                  method = "pearson")$estimate,
                                     simplify = T),
                  "spearman" = sapply(rownames(predicted),
                                      function(gene_name) cor.test(x = centered_data[gene_name,],
                                                                   y = predicted[gene_name,],
                                                                   method = "spearman")$estimate,
                                      simplify = T))
    
    names(cor) <- sapply(names(cor), function(x) strsplit(x,"\\.")[[1]][1])
    
    if(COR == "pearson"){
      saveRDS(cor, paste0(PRED_DIR,"/",gsub("net_predictions",
                                            "correlations",f[i])))      
    } else if(COR == "spearman"){
      saveRDS(cor, paste0(PRED_DIR,"/",gsub("net_predictions",
                                            "correlations_spearman",f[i])))
    }
  }
}



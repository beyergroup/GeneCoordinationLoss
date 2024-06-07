# Compute correlations between network predictions and original data

args <- commandArgs(trailingOnly = T)
PATTERN = args[1]
EXPR_DIR = args[2]
PRED_DIR = args[3]
MODE = args[4]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")

files <- list.files(PRED_DIR, pattern = PATTERN)
files <- files[grep("sampled_centered_net_predictions",files)]

data_files <- list.files(EXPR_DIR, pattern = PATTERN)
data_files <- data_files[grep("sampled_centered_data",data_files)]


for(i in 1:length(files)){
  
  predicted <- as.matrix(ReadRDS(paste0(PRED_DIR,"/",files[i])))
  centered_data <- as.matrix(ReadRDS(paste0(EXPR_DIR,"/",data_files[i])))
  centered_data <- DataENSGToSymbol(centered_data)
  centered_data <- centered_data[rownames(predicted),]
  
  if(sum(is.na(predicted)) > 0)
    break
  
  if(MODE == "Pearson"){
    
    message("Computing Pearson correlations")
    
    cor <- sapply(rownames(predicted),
                  function(gene_name) cor.test(x = centered_data[gene_name,],
                                               y = predicted[gene_name,],
                                               method = "pearson")$estimate,
                  simplify = T)
    names(cor) <- sapply(names(cor), function(x) strsplit(x,"\\.")[[1]][1])
    
    WriteRDS(cor, paste0(PRED_DIR,"/",gsub("net_predictions",
                                           "correlations", files[i])))
    rm(cor); gc()
    
  } else if(MODE == "Spearman"){
    
    message("Computing Spearman correlations")
    
    cor <- sapply(rownames(predicted),
                  function(gene_name) cor.test(x = centered_data[gene_name,],
                                               y = predicted[gene_name,],
                                               method = "spearman")$estimate,
                  simplify = T)
    names(cor) <- sapply(names(cor), function(x) strsplit(x,"\\.")[[1]][1])
    
    WriteRDS(cor, paste0(PRED_DIR,"/",gsub("net_predictions",
                                           "correlations_spearman",files[i])))
    rm(cor); gc()
  }
}
rm(centered_data,predicted); gc()

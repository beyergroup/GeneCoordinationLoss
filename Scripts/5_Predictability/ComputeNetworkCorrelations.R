# Compute correlations between network predictions and original data

args <- commandArgs(trailingOnly = T)
EXPR_DIR = args[1]
PRED_DIR = args[2]

source("Scripts/functions.R")

files <- list.files(PRED_DIR, pattern = "centered_net_predictions.rds")
files <- files[grep("young|old",files)]
# files <- files[grep("v7",files)]
data_files <- list.files(EXPR_DIR, pattern = "centered_data.rds")
# data_files <- list.files(PRED_DIR, pattern = "v7_sampled_centered_data.rds")
data_files <- data_files[-grep("v7",data_files)]

for(i in 1:length(files)){
  
  predicted <- as.matrix(readRDS(paste0(PRED_DIR,"/",files[i])))
  centered_data <- as.matrix(readRDS(paste0(EXPR_DIR,"/",data_files[i])))
  centered_data <- DataENSGToSymbol(centered_data)
  centered_data <- centered_data[rownames(predicted),]
  
  if(sum(is.na(predicted)) > 0)
    break
  
  cor <- sapply(rownames(predicted), function(gene_name) cor.test(x = centered_data[gene_name,],
                                                                  y = predicted[gene_name,],
                                                                  method = "pearson")$estimate, simplify = T)
  names(cor) <- sapply(names(cor), function(x) strsplit(x,"\\.")[[1]][1])
  
  saveRDS(cor, paste0(PRED_DIR,"/",gsub("net_predictions","correlations",files[i])))
  rm(cor)
  
  cor <- sapply(rownames(predicted), function(gene_name) cor.test(x = centered_data[gene_name,],
                                                                  y = predicted[gene_name,],
                                                                  method = "spearman")$estimate, simplify = T)
  names(cor) <- sapply(names(cor), function(x) strsplit(x,"\\.")[[1]][1])
  
  saveRDS(cor, paste0(PRED_DIR,"/",gsub("net_predictions","correlations_spearman",files[i])))
}

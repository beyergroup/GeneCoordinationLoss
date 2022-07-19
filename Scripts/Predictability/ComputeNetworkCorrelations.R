# Compute correlations between network predictions and original data

args <- commandArgs(trailingOnly = T)

NET_FILE = "Data/Networks/Human/stabsel_network_Hs_filtered.rds"
EXPR_DIR = "Outputs/Tabula_Sapiens/SampledSCTNorm/"
PRED_DIR = "Outputs/Human_Network/stabsel/Predictability/Tabula_Sapiens/"

setwd("../../")

source("Scripts/functions.R")

files <- list.files(PRED_DIR, pattern = "net_predictions.rds")
data_files <- list.files(EXPR_DIR, pattern = "centered_data.rds")

for(i in 1:length(files)){
  
  predicted <- as.matrix(readRDS(paste0(PRED_DIR,files[i])))
  centered_data <- as.matrix(readRDS(paste0(EXPR_DIR,data_files[i])))
  
  centered_data <- centered_data[rownames(predicted),]
  
  if(sum(is.na(predicted)) > 0)
    break
  
  cor <- sapply(rownames(predicted), function(gene_name) cor.test(x = centered_data[gene_name,],
                                                                  y = predicted[gene_name,],
                                                                  method = "pearson")$estimate, simplify = T)
  names(cor) <- sapply(names(cor), function(x) strsplit(x,"\\.")[[1]][1])
  
  saveRDS(cor, paste0(PRED_DIR,gsub("net_predictions","correlations",files[i])))
}

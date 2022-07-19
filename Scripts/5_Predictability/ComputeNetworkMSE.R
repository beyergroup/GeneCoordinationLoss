# Compute MSE per age group and tissue

EXPR_DIR = "GTEx_Networks/AgeTissue_Networks/Outputs"
PRED_DIR = "Outputs/5_Predictability"

source("Scripts/functions.R")

files <- list.files(PRED_DIR, pattern = "net_predictions.rds")
files <- files[grep("young|old",files)]
data_files <- list.files(EXPR_DIR, pattern = "centered_data.rds")

for(i in 1:length(files)){
  
  predicted <- as.matrix(readRDS(paste0(PRED_DIR,"/",files[i])))
  centered_data <- as.matrix(readRDS(paste0(EXPR_DIR,"/",data_files[i])))
  centered_data <- DataENSGToSymbol(centered_data)
  centered_data <- centered_data[rownames(predicted),]
  
  if(sum(is.na(predicted)) > 0)
    break
  
  mse <- sapply(rownames(predicted),
                function(gene_name) ComputeNMSE(pred = predicted[gene_name,],
                                                obs = centered_data[gene_name,]))
  
  saveRDS(mse, paste0(PRED_DIR,"/",gsub("net_predictions","NMSE",files[i])))
  
  rm(predicted, centered_data, mse); gc()
}

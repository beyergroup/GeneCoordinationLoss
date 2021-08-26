# Compute delta matrices

args = commandArgs(trailingOnly=TRUE)
NET = args[1]

PREDICTION_FOLDER = paste0("Outputs/Human_Network/",NET,"/Predictability/AgeTissue")
DATA_FOLDER = "GTEx_Networks/AgeTissue_Networks/Outputs"

setwd("../")

files <- list.files(PREDICTION_FOLDER, pattern = "net_predictions.rds", full.names = T)
data_files <- list.files(DATA_FOLDER, pattern = "_sampled_centered_data.rds", full.names = T)

for(i in 1:length(files)){
  
  predicted <- readRDS(files[i])
  centered_data <- readRDS(data_files[i])
  
  # convert rownames to gene symbols
  conversion_table <- read.delim("Resources/ensembl_idversion_GTExDESeq2_symbolChrStart.txt")
  rownames(centered_data) <- sapply(rownames(centered_data),
                                    function(c) strsplit(c, split = "\\.")[[1]][1])
  rnames <- conversion_table[match(rownames(centered_data), conversion_table$ensembl_gene_id),"symbol"]
  centered_data <- centered_data[!is.na(rnames),]
  rownames(centered_data) <- rnames[!is.na(rnames)]
  rm(conversion_table,rnames); gc()
  
  centered_data <- centered_data[rownames(predicted),]
  
  if(sum(is.na(predicted)) > 0)
    break
  
  delta <- abs(predicted-centered_data)
  
  saveRDS(delta, gsub("net_predictions","delta",files[i]))
}

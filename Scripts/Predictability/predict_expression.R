# Prediction of expression based on the networks

# NET_FILE = "Data/Networks/Human/stabsel_pcclasso_network_Hs.rds"
# PATH = "Human_Network/stabsel_pcclasso/"
# NET_FILE = "Data/Networks/Human/stabsel_pcclasso_network_Hs_filtered.rds"
# PATH = "Human_Network/stabsel_pcclasso_filter01/"
# NET_FILE = "Data/Networks/Human/stabsel_network_Hs.rds"
# PATH = "Human_Network/stabsel/"
NET_FILE = "Data/Networks/Human/stabsel_randomized_network_Hs.rds"
PATH = "Human_Network/stabsel_randomized/"

DATA_FOLDER = "GTEx_Networks/Tissue_Networks/Outputs"
# DATA_FOLDER = "GTEx_Networks/Age_Networks/Outputs"
# DATA_FOLDER = "GTEx_Networks/AgeTissue_Networks/Outputs"

source("functions.R")

setwd("../")

net <- as.matrix(readRDS(NET_FILE))
data_files <- list.files(DATA_FOLDER, pattern = "_sampled_centered_data.rds", full.names = T)

for(file in data_files){
  
  centered_data <- readRDS(file)
  
  prefix <- tail(strsplit(file, "/")[[1]],1)
  
  # convert rownames to gene symbols
  conversion_table <- read.delim("Resources/ensembl_idversion_GTExDESeq2_symbolChrStart.txt")
  rownames(centered_data) <- sapply(rownames(centered_data),
    function(c) strsplit(c, split = "\\.")[[1]][1])
  rnames <- conversion_table[match(rownames(centered_data), conversion_table$ensembl_gene_id),"symbol"]
  centered_data <- centered_data[!is.na(rnames),]
  rownames(centered_data) <- rnames[!is.na(rnames)]
  rm(conversion_table,rnames); gc()
  
  pred <- PredictNet(net, centered_data, maxiter = 1)
  saveRDS(pred, paste0("Outputs/",PATH,"Predictability/Tissue/",
                       gsub("sampled_centered_data","sampled_net_predictions",prefix)))
  
  rm(centered_data,pred,prefix); gc()
}

rm(net); gc()

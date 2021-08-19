# Prediction of expression based on the networks

args <- commandArgs(trailingOnly = T)
NET = args[1]
TYPE = args[2]
RAND = as.logical(args[3])

NET = "stabsel" # "stabsel_pcclasso"
TYPE = "Tissue" # "AgeTissue" # "Age"
RAND = F # T

if(RAND){
  net_file = paste0("Data/Networks/Human/",NET,"_randomized_network_Hs_filtered.rds")
  path = paste0("Human_Network/",NET,"_randomized/","Predictability/",TYPE,"/")
} else{
  net_file = paste0("Data/Networks/Human/",NET,"_network_Hs_filtered.rds")
  path = paste0("Human_Network/",NET,"/Predictability/",TYPE,"/")
}
dir.create(paste0("Outputs/",path))

data_folder = paste0("GTEx_Networks/",TYPE,"_Networks/Outputs")

source("Scripts/functions.R")


net <- as.matrix(readRDS(net_file))
data_files <- list.files(data_folder, pattern = "_sampled_centered_data.rds", full.names = T)

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
  saveRDS(pred, paste0("Outputs/",path,
                       gsub("sampled_centered_data","sampled_net_predictions",prefix)))
  
  rm(centered_data,pred,prefix); gc()
}

rm(net); gc()

# Predict target gene expression

args = commandArgs(trailingOnly = T)
NET = args[1]
PATTERN = args[2]
EXPR_DIR = args[3]
OUTDIR = args[4]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")

dir.create(OUTDIR)

net <- as.matrix(ReadRDS(paste0("Outputs/0_Preprocessing/",NET,"_network_Hs.rds")))

if(all(grepl("ENSG",rownames(net)))){
  rownames(net) <- VectorENSGToSymbol(rownames(net))
  colnames(net) <- VectorENSGToSymbol(colnames(net))
}

data_files <- list.files(EXPR_DIR, pattern = PATTERN)
data_files <- data_files[grep("sampled_centered_data",data_files)]

for(file in data_files){
  
  centered_data <- ReadRDS(paste0(EXPR_DIR,"/",file))
  centered_data <- DataENSGToSymbol(centered_data)
  
  pred <- PredictNet(net, as.matrix(centered_data), maxiter = 1)
  WriteRDS(pred, paste0(OUTDIR,"/",
                        gsub("_data","_net_predictions",file)))
  
  rm(centered_data,pred); gc()
}

rm(net,data_files); gc()

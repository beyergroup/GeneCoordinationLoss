# Get gene metrics in each tissue subset

WDIR = "/data/public/adesous1/GeneCorrelation/"
setwd(WDIR)

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")

files <- list.files("GTEx_Networks/Tissue_Networks/Outputs",
                    pattern = "sampled_data.rds", full.names = T)


for(file in files){
  
  gtex <- readRDS(file)
  
  # convert to gene symbol
  gtex <- DataENSGToSymbol(gtex, remove_dup = T)
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  # "logFC" - residuals after batch correction
  LFCs <- rowMeans(gtex)
  saveRDS(LFCs, paste0("Outputs/GTEx/logFC/",tissue,
                       "_sampled_meanLFC.rds"))
  rm(LFCs); gc()
  
  
  # variance
  var <- apply(gtex, 1, var)
  saveRDS(var, paste0("Outputs/GTEx/logFC/",tissue,
                       "_sampled_var.rds"))
  rm(var,gtex); gc()
}

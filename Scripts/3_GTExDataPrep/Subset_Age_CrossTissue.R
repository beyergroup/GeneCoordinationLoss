# Subset GTEx data

args <- commandArgs(trailingOnly = TRUE)
DATA_DIR = args[1]
SUBSET_DIR = args[2]

library(ggplot2)
source("Scripts/functions.R")
source("Scripts/3_GTExDataPrep/params.R")

metadata <- ReadRDS(paste0("Outputs/3_GTExDataPrep/metadata_restrictedaccess.rds"))


for(age_group in c("20-29","30-39","40-49","50-59","60-69","70-79")){
  
  data <- c()
  
  for(tissue in MAIN_TISSUES){

    if(!file.exists(paste0(SUBSET_DIR,"/",tissue,"_",
                           age_group,"_sampled_data.rds")))
      next
    
    tmp <- ReadRDS(paste0(SUBSET_DIR,"/",tissue,"_",
                          age_group,"_sampled_data.rds"))
    
    data <- cbind(data, tmp[,sample(1:ncol(tmp),2)])
  }
  
  WriteRDS(data, paste0("Outputs/3_GTExDataPrep/Subset_Data/CrossTissue_",
                        age_group,"_sampled_data.rds"))
}


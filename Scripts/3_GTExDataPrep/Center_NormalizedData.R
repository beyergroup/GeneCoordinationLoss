args <- commandArgs(trailingOnly = TRUE)
WDIR = args[1]

source("Scripts/functions.R")

files <- list.files(WDIR, pattern = "sampled_data.rds")

for(file in files){
  
  data <- ReadRDS(paste0(WDIR,"/",file))
  
  # method <- gsub(".rds","",tail(strsplit(file,"_")[[1]],1))
  
  # center
  centers <- rowMeans(data)
  data <-  sweep(data, 1, centers, "-")
  
  # save
  WriteRDS(data, paste0(WDIR,"/",gsub("sampled","sampled_centered",file)))
}

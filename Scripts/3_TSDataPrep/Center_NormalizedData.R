
WDIR = "3_TSDataPrep"

source("Scripts/functions.R")
source("Scripts/3_TSDataPrep/params.R")

files <- list.files(paste0("Outputs/",WDIR,"/Normalized/Subset/"),
                    pattern = ".rds")
files <- files[-grep("quantile",files)]

for(file in files){
  
  data <- ReadRDS(paste0("Outputs/",WDIR,"/Normalized/Subset/",file))
  
  method <- gsub(".rds","",tail(strsplit(file,"_")[[1]],1))
  
  # remove genes not quantified in enough cells
  data[is.na(data)] <- 0
  data <- data[rowSums(data != 0) >= QCELL_THRE[method],]
  message(nrow(data)," genes included")
  
  # center
  centers <- rowMeans(data)
  data <-  sweep(data, 1, centers, "-")
  
  # save
  WriteRDS(data, paste0("Outputs/",WDIR,"/Centered/",
                        gsub(".rds","_centered.rds",
                             gsub("alldetected","quantified",file))))
}

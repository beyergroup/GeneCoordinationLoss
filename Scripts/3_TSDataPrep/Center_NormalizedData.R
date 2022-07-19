
WDIR = "3_TSDataPrep"

source("Scripts/functions.R")
source("Scripts/3_TSDataPrep/params.R")

files <- list.files(paste0("Outputs/",WDIR,"/Normalized"),
                    pattern = ".rds")
files <- files[-grep("quantile",files)]

for(file in files){
  
  sce <- ReadRDS(paste0("Outputs/",WDIR,"/Normalized/",file))
  
  method <- gsub(".rds","",tail(strsplit(file,"_")[[1]],1))
  
  # check if cell number is high enough
  if(ncol(sce) < NCELL_THRE[method]){
    next
  }
  
  # remove genes not quantified in enough cells
  assay(sce,"SCT")[is.na(assay(sce,"SCT"))] <- 0
  sce <- sce[rowSums(assay(sce,"SCT") != 0) >= QCELL_THRE[method],]
  message(nrow(sce)," genes included")
  
  # center
  centers <- rowMeans(assay(sce,"SCT"))
  data <-  sweep(assay(sce,"SCT"), 1, centers, "-")
  
  # save
  WriteRDS(data, paste0("Outputs/",WDIR,"/Centered/",
                        gsub(".rds","_centered.rds",
                             gsub("alldetected","quantified",file))))
}

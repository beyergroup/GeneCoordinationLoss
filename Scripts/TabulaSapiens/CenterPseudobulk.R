# Center pseudobulk

DIR = "../../Outputs/Tabula_Sapiens/Pseudobulk/"

files <- list.files(DIR, pattern = "list.rds")

for(file in files){
  
  list <- readRDS(paste0(DIR,file))
  
  for(type in c("Sum","Mean","Median")){
    dir.create(paste0(DIR,type))
    means <- rowMeans(list[[type]])
    centered <- sweep(list[[type]], 1, means, "-")
    saveRDS(centered, paste0(DIR,type,"/",gsub("_list","",file)))
  }
}

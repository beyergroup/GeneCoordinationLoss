# Correlation all genes vs genes

source("Scripts/functions.R")

# GTEx -------------------------------------------------------------------------

gtex_files <- list.files("GTEx_Networks/Tissue_Networks/Outputs",
                         pattern = "sampled_centered_data.rds")

for(file in gtex_files){
  
  data <- ReadRDS(paste0("GTEx_Networks/Tissue_Networks/Outputs/",file))
  data <- DataENSGToSymbol(data, remove_dup = T)
  
  cor <- cor(t(data))
  
  WriteRDS(cor, paste0("Outputs/0.5_Correlation/",
                       gsub("sampled_centered_data","correlation",file)))
  
  rm(data,cor); gc()
}
rm(file,gtex_files); gc()


# Tabula Sapiens ---------------------------------------------------------------

sc_files <- list.files("Outputs/3_TSDataPrep/Normalized/Subset",
                       pattern = "quantile")

for(file in sc_files){
  
  data <- ReadRDS(paste0("Outputs/3_TSDataPrep/Normalized/Subset/",file))
  cor <- cor(t(data), use = "pairwise.complete.obs")
  
  WriteRDS(cor, paste0("Outputs/0.5_Correlation/",
                       gsub("quantile","correlation",file)))
  
  rm(data,cor); gc()
}
rm(file,sc_files); gc()
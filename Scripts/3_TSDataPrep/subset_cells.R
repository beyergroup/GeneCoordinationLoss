
source("Scripts/functions.R")

METHOD = "smartseq2"
files <- list.files("Outputs/3_TSDataPrep/Normalized/Subset/",
                    pattern = paste0(METHOD,".rds"))

subset_cells <- list()

for(file in files){
  cell_type <- paste(strsplit(file,"_")[[1]][-c(1,2,length(strsplit(file, 
                                                                    "_")[[1]]))],
                     collapse = "_")
  sce <- ReadRDS(paste0("Outputs/3_TSDataPrep/Normalized/Subset/",file))
  subset_cells[[cell_type]] <- colnames(sce)
}

WriteRDS(subset_cells, paste0("Outputs/3_TSDataPrep/subset_cells_",
                              METHOD,".rds"))

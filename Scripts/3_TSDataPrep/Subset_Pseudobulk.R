# Subset pseudobulk per cell-type/tissue combination

args <- commandArgs(trailingOnly = T)
WDIR = args[1]
METHOD = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")
source(paste0("Scripts/",WDIR,"/params.R"))


norm <- ReadRDS(paste0("Outputs/3_TSDataPrep/",METHOD,
                       "/Pseudobulk/pseudobulk_DESeq2norm.rds"))
metadata <- ReadRDS(paste0("Outputs/3_TSDataPrep/",METHOD,
                           "/Pseudobulk/pseudobulk_metadata.rds"))

# combinations <- unique(paste(metadata$Tissue, metadata$CellType, sep = "-"))

for(tissue in unique(metadata$Tissue)){
  
  tmp <- subset(metadata, Tissue == tissue)
  
  for(cell_type in unique(tmp$CellType)){
    
    if(sum(tmp$CellType == cell_type) >= 5){
      
      # subset to selection
      samples <- subset(tmp, CellType == cell_type)$Sample
      data <- norm[,sample(samples)[1:5]]
      WriteRDS(data, paste0("Outputs/",WDIR,"/",METHOD,"/Pseudobulk/Subsets/",
                            tissue,"_",cell_type,".rds"))
    }
  }
}

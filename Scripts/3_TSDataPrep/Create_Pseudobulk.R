# Create pseudobulk from sum of reads across cells of same cell type, organ and donor

args <- commandArgs(trailingOnly = T)
WDIR = args[1]
METHOD = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(SingleCellExperiment)
source("Scripts/functions.R")
source(paste0("Scripts/",WDIR,"/params.R"))


# files <- list.files(paste0("Outputs/",WDIR,"/",METHOD))
files <- list.files(paste0("Outputs/",WDIR,"/Filters"),
                    pattern = "alldetected_genes", full.names = TRUE)
files <- files[grep(METHOD,files)]

pseudobulk <- c()
metadata <- data.frame()

for(file in files){
  
  sce <- ReadRDS(file)
  
  donors <- as.character(unique(sce@colData$donor))
  # cell_types <- as.character(unique(sce@colData$cell_ontology_class))
  cell_type <- as.character(unique(sce@colData$cell_ontology_class))
  tissue <- as.character(unique(sce@colData$organ_tissue))
  
  for(donor in donors){
    
    tmp <- sce[,(sce@colData$donor == donor) &
                 (sce@colData$cell_ontology_class == cell_type)]
    
    if(ncol(tmp) < NCELL_THRE[METHOD]){
      next
    }
    
    # subsample
    # tmp <- tmp[,sample(1:ncol(tmp), NCELL_THRE[METHOD])]
    tmp <- tmp[,sample(colnames(tmp))]
    
    indexes <- seq(from = 1, by = PB_CELLN[METHOD],
                   length.out = floor(ncol(tmp)/PB_CELLN[METHOD]))
    
    for(i in indexes){
      # sum over cells for pseudobulk
      pseudobulk <- cbind(pseudobulk,
                          rowSums(assay(tmp,"decontXcounts")[,i:(i+PB_CELLN[METHOD]-1)]))
      metadata <- rbind.data.frame(metadata,
                                   data.frame("Tissue" = tissue,
                                              "Donor" = donor,
                                              "CellType" = cell_type))
    }
  }
  rm(sce,tmp); gc()
}

colnames(pseudobulk) <- paste0("PB",1:ncol(pseudobulk))
metadata$Sample <- paste0("PB",1:ncol(pseudobulk))
metadata$Method <- METHOD

WriteRDS(pseudobulk, paste0("Outputs/",WDIR,"/",METHOD,
                            "/Pseudobulk/pseudobulk_data.rds"))
WriteRDS(metadata, paste0("Outputs/",WDIR,"/",METHOD,
                          "/Pseudobulk/pseudobulk_metadata.rds"))

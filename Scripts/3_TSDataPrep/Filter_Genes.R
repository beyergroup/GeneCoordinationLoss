# Gene filters

args <- commandArgs(trailingOnly = T)
WDIR = args[1]
METHOD = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")
source(paste0("Scripts/",WDIR,"/params.R"))
library(SingleCellExperiment)
library(ggplot2)
library(cowplot)
library(tidyverse)


# Remove genes not detected in any cells of cell type
files <- list.files(paste0("Outputs/",WDIR,"/Filters"),
                    recursive = T, full.names = T, pattern = METHOD)
files <- files[grep("cellfiltered_",files)]

for(file in files){
  
  tissue <- gsub(tail(strsplit(file,"_")[[1]],1), pattern = ".rds", replacement = "")
  
  sce <- ReadRDS(file)
  
  # separate cell types
  cell_types <- as.character(unique(sce@colData$cell_ontology_class))
  
  for(cell_type in cell_types){
    
    # subset
    tmp_sce <- sce[,sce@colData$cell_ontology_class == cell_type]
    
    # remove undetected genes
    tmp_sce <- tmp_sce[rowSums(assay(tmp_sce,"decontXcounts") != 0) != 0,]
    
    # save
    ct <- gsub(cell_type, pattern = " ", replacement = "_")
    saveRDS(tmp_sce, paste0("Outputs/",WDIR,"/Filters/alldetected_genes_",
                            tissue,"_",ct,"_",METHOD,".rds"), version = 2)
    rm(tmp_sce,ct); gc()
  }
  
  rm(sce,cell_types,tissue); gc()
}

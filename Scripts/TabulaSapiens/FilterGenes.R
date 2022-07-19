# Gene filters

.libPaths("../../Resources/Rlibs/R-4.0.3/")

# Remove genes not detected in any cells of cell type
files <- list.files("../../Outputs/Tabula_Sapiens/CellFiltered",
                    recursive = T, full.names = T)
for(file in files){
  
  tissue <- gsub(tail(strsplit(file,"_")[[1]],1), pattern = ".rds", replacement = "")
  method <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  sce <- readRDS(file)
  
  # separate cell types
  cell_types <- as.character(unique(sce@colData$cell_ontology_class))
  
  for(cell_type in cell_types){
    
    # subset
    tmp_sce <- sce[,sce@colData$cell_ontology_class == cell_type]
    
    # remove undetected genes
    tmp_sce <- tmp_sce[rowSums(assay(tmp_sce,"decontXcounts") != 0) != 0,]
    
    # save
    ct <- gsub(cell_type, pattern = " ", replacement = "_")
    saveRDS(tmp_sce, paste0("../../Outputs/Tabula_Sapiens/GeneFiltered/all_detected_genes/",
                            tissue,"_",ct,"_",method,".rds"), version = 2)
    rm(tmp_sce,ct); gc()
  }
  
  rm(sce,cell_types,tissue,method); gc()
}

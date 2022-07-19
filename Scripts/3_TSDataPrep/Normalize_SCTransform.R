# sctransform normalization

args <- commandArgs(trailingOnly = T)
WDIR = args[1]
METHOD = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")
library(Seurat)
library(glmGamPoi)


files <- list.files(paste0("Outputs/",WDIR,"/Filters"),
                    pattern = METHOD, full.names = T)
files <- files[grep("alldetected_genes",files)]

for(file in files){
  
  sce <- ReadRDS(file)
  
  # convert to Seurat
  seu <- CreateSeuratObject(counts = as.matrix(assay(sce,"decontXcounts")),
                            meta.data = as.data.frame(sce@colData))
  
  # Apply SCTransform
  seu <- SCTransform(seu, method = "glmGamPoi", vst.flavor = "v2")
  
  # Convert back to SingleCellExperiment
  sce <- sce[rownames(seu),]
  assay(sce,"SCT") <- as(seu@assays$SCT@data, "dgTMatrix")
  
  saveRDS(sce, gsub("/Filters/","/Normalized/",file))
  rm(sce,seu); gc()
}

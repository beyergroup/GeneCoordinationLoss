# sctransform normalization

.libPaths("../../Resources/Rlibs/R-4.0.3/")
library(Seurat)
library(glmGamPoi)

INDIR = "../../Outputs/Tabula_Sapiens/GeneFiltered/all_detected_genes"
OUTDIR = "../../Outputs/Tabula_Sapiens/SCTransformNormalized/all_detected_genes"
cluster = F
# PATTERN = "Lung"
PATTERN = NULL

files <- list.files(INDIR, pattern = PATTERN)


for(file in files[263:319]){
  
  sce <- readRDS(paste0(INDIR,"/",file))
  
  if(ncol(sce) < 100){
    next
  }
  
  # convert to Seurat
  seu <- CreateSeuratObject(counts = as.matrix(assay(sce,"decontXcounts")),
                            meta.data = as.data.frame(sce@colData))
  
  # Apply SCTransform
  seu <- SCTransform(seu, method = "glmGamPoi")
  
  # Convert back to SingleCellExperiment
  sce <- sce[rownames(seu),]
  assay(sce,"SCT") <- as(seu@assays$SCT@data, "dgTMatrix")
  
  saveRDS(sce, paste0(OUTDIR,"/",file))
  rm(sce,seu); gc()
}

# scran normalization

.libPaths("../../Resources/Rlibs/R-4.0.3/")
library(scran)
library(scuttle)


INDIR = "../../Outputs/Tabula_Sapiens/GeneFiltered/all_detected_genes"
OUTDIR = "../../Outputs/Tabula_Sapiens/scranNormalized/all_detected_genes"
cluster = F
PATTERN = "Pancreas"

files <- list.files(INDIR, pattern = PATTERN)

for(file in files){
  
  sce <- readRDS(paste0(INDIR,"/",file))
  sce <- computeSumFactors(sce, assay.type = "decontXcounts")
  sce <- logNormCounts(sce, assay.type = "decontXcounts")
  
  saveRDS(sce, paste0(OUTDIR,"/",file))
  rm(sce); gc()
}

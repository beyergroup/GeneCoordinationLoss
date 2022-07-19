# Split tissue data by sequencing method and tissue

args <- commandArgs(trailingOnly = T)
INDIR = args[1]
WDIR = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(zellkonverter)
library(SingleCellExperiment)

files <- list.files(INDIR, full.names = TRUE, pattern = ".h5ad")
files <- files[-grep(".zip",files)]
files <- files[-grep("endothelial|epithelial|immune|stromal",files)]

dir.create(paste0("Outputs/",WDIR,"/10X"))
dir.create(paste0("Outputs/",WDIR,"/smartseq2"))

for(file in files){
  
  sce <- readH5AD(file)
  
  # separate by sequencing method
  for(m in unique(sce@colData$method)){
    saveRDS(sce[, sce@colData$method == m],
            paste0("Outputs/",WDIR,"/",m,"/",
                   unique(sce@colData$organ_tissue),
                   ".rds"))
  }
}

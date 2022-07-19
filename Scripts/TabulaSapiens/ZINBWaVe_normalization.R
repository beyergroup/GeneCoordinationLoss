# ZINB-WaVe normalization

args <- commandArgs(trailingOnly = T)
INDIR = args[1]
OUTDIR = args[2]
cluster = as.logical(args[3])
PATTERN = args[4]

INDIR = "../../Outputs/Tabula_Sapiens/GeneFiltered/all_detected_genes"
OUTDIR = "../../Outputs/Tabula_Sapiens/ZBNormalized/all_detected_genes"
cluster = F
PATTERN = "Vasculature"

if(cluster){
  message("Using cluster version (R 3.4.4)")
  # library(zinbwave, lib.loc = "/cellnet/Ronja/R/")
  .libPaths(new = "../../Resources/Rlibs/R-3.4.4/")
  devtools::load_all("../../Resources/Rlibs/R-3.4.4/pcaPP/")
  devtools::load_all("../../Resources/Rlibs/R-3.4.4/pspline/")
  devtools::load_all("../../Resources/Rlibs/R-3.4.4/edgeR/")
  devtools::load_all("../../Resources/Rlibs/R-3.4.4/genefilter/")
  devtools::load_all("../../Resources/Rlibs/R-3.4.4/softImpute/")
  devtools::load_all("../../Resources/Rlibs/R-3.4.4/zinbwave/")
} else{
  .libPaths(new = "../../Resources/Rlibs/R-4.0.3/")
  library(zinbwave)
}

# ZINBWaVe observational weights

files <- list.files(INDIR, pattern = PATTERN)

ZINBWaVeCluster <- function(sce, k = 2, cluster){
  
  library(SummarizedExperiment)
  
  if(cluster){
    options(matprod="internal")
    BiocParallel::register(BiocParallel::SerialParam())
  }
  
  assay(sce,"X") <- NULL
  assay(sce,"raw_counts") <- NULL
  
  assay(sce,"decontXcounts") <- as(assay(sce,"decontXcounts"), "dgCMatrix")
  fit <- zinbFit(Y = sce,
                 X = "~ n_counts_UMIs",
                 K = k,
                 verbose = T,
                 which_assay = "decontXcounts")
  
  assay(sce,"decontXcounts") <- as(assay(sce,"decontXcounts"), "matrix")
  z <- zinbwave(sce,
                fitted_model = fit,
                verbose = T,
                observationalWeights = TRUE,
                normalizedValues = TRUE,
                residuals = TRUE,
                which_assay = "decontXcounts")
  
  return(z)
}

k = 2


for(file in files){
  
  cat(file,"\n")
  data <- readRDS(paste0(INDIR,"/",file))
  data <- ZINBWaVeCluster(data, k = k, cluster = cluster)
  saveRDS(data,
          # file = paste0(OUTDIR,"/",file))
          file = gsub(paste0(OUTDIR,"/",file), pattern = ".rds", replacement = "_LibSize.rds"))
  rm(data); gc()
}

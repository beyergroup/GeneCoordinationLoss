# Evaluate normalization stragies

NAME = "SCTransform"
# ASSAY = "logcounts"
ASSAY = "SCT"
TISSUES = c("Lung","Pancreas","Vasculature")

.libPaths("../../Resources/Rlibs/R-4.0.3/")
library(ggpubr)
library(Seurat)
library(SingleCellExperiment)
source("functions.R")

tissue = "Lung"
INDIR = paste0("../../Outputs/Tabula_Sapiens/",NAME,
               "Normalized/all_detected_genes")

files <- list.files(INDIR, pattern = tissue)

for(file in files){
  
  # Extract cell type and method
  method <- gsub(paste(strsplit(file,"_")[[1]][grep("10X|smartseq2",
                                                    strsplit(file,"_")[[1]]):length(strsplit(file,"_")[[1]])],
                       collapse = "_"),
                 pattern = ".rds", replacement = "")
  cell_type <- paste(strsplit(file,"_")[[1]][2:(grep("10X|smartseq2",strsplit(file,"_")[[1]])-1)],
                     collapse = " ")
  
  sce <- readRDS(paste0(INDIR,"/",file))
  
  plots <- list()
  # Expression profiles of random cells
  plots[[1]] <- NormEval_ExprScatter(sce, assay = ASSAY,
                       title = c(paste0(tissue,", ",cell_type),
                                 paste0("- ",method)))
  
  # Trend with LibSize by quantile
  plots[[2]] <- NormEval_LibTrend(sce, assay = ASSAY,
                    raw_assay = "decontXcounts", Q = 5,
                    title = c(paste0(tissue,", ",cell_type),
                              paste0("- ",method)))
  
  png(paste0("../../Plots/Tabula_Sapiens/QC/",tissue,"_",
             gsub(cell_type,pattern = " ",replacement = "_"),
             method,"_",NAME,"_corr.png"), width = 700, height = 500)
  print(ggarrange(plotlist = plots, widths = c(3,1.9), align = "h"))
  dev.off()
  
  rm(sce,method,cell_type,plots); gc()
}

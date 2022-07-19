# DESeq2 normalization

args <- commandArgs(trailingOnly = T)
WDIR = args[1]
METHOD = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(SingleCellExperiment)
library(DESeq2)
source("Scripts/functions.R")
source(paste0("Scripts/",WDIR,"/params.R"))


pseudobulk <- ReadRDS(paste0("Outputs/3_TSDataPrep/",METHOD,
                             "/Pseudobulk/pseudobulk_data.rds"))
metadata <- ReadRDS(paste0("Outputs/3_TSDataPrep/",METHOD,
                           "/Pseudobulk/pseudobulk_metadata.rds"))


dds <- DESeq2::DESeqDataSetFromMatrix(countData = pseudobulk,
                                      colData = metadata, design = ~ 1)
dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
norm <- counts(dds, normalized = T)
WriteRDS(norm, paste0("Outputs/",WDIR,"/",METHOD,
                      "/Pseudobulk/pseudobulk_DESeq2norm.rds"))

# type = "poscounts" because the default way to compute the geometric means does
# not allow for any 0 in the matrix
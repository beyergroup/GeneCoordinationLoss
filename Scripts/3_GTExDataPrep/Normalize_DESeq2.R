
args <- commandArgs(trailingOnly = TRUE)
OUT_DIR = args[1]

library(DESeq2, include.only = c("counts"))

## 01 - DESeq2 normalization

annotation <- readRDS(paste0(OUT_DIR,"/annotation.rds"))
gtex <- readRDS(paste0(OUT_DIR,"/raw_gtex.rds"))

gtex <- DESeqDataSetFromMatrix(countData = gtex, colData = annotation, design = ~ 1)
gtex <- estimateSizeFactors(gtex)

gtex <- counts(gtex, normalized = T)
SaveRDS(gtex, paste0(OUT_DIR,"/DESeq2_gtex.rds"))
SaveRDS(annotation, paste0(OUT_DIR,"/DESeq2_annotation.rds"))

rm(gtex,annotation); gc()


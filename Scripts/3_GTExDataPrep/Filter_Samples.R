# Filter out low quality and cell line samples from GTEx

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR = args[1]

gtex <- readRDS(paste0(OUT_DIR,"/DESeq2_gtex.rds"))
annotation <- readRDS(paste0(OUT_DIR,"/DESeq2_annotation.rds"))


# Filter Samples ---------------------------------------------------------------------
# RIN and cells
annotation <- annotation[(annotation$SMRIN >= 6) & (!grepl("Cells",annotation$SMTSD)),]
gtex <- gtex[,annotation$SAMPID]
SaveRDS(annotation, paste0(OUT_DIR,"/sampleFiltered_annotation.rds"))
SaveRDS(gtex, paste0(OUT_DIR,"/sampleFiltered_gtex.rds"))
#I filtered samples before than genes so I would not choose genes based on samples that I would further remove.
#-------------------------------------------------------------------------------------

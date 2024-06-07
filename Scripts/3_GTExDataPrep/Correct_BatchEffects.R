# Correct batch effects per tissue (SMTSD)

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR = args[1]

annotation <- readRDS(paste0(DATA_DIR,"sampleFiltered_annotation.rds"))
gtex <- readRDS(paste0(DATA_DIR,"sample_and_gene_biotype_and_mean_filtered_gtex.rds")) #rows=genes cols=samples

unique(annotation$SMTSD)
unique(annotation$SMTS)

gtex_final <- c()
annotation_final <- c()
ltt <- unique(annotation$SMTSD) 
for (tt in ltt) {
  gtex_sub <- gtex[,annotation$SAMPID[annotation$SMTSD==tt]]
  
  # reduce annotation to known covariates of interest: ischemic time, experimental batch and death type
  annotation_sub <- annotation[annotation$SMTSD==tt, c("SAMPID","SUBJID","SMTSISCH","SMGEBTCH","DTHHRDY", "AGE","SMTSD","SEX")]
  annotation_sub$DTHHRDY <- as.factor(annotation_sub$DTHHRDY)
  annotation_sub$SEX <- as.factor(annotation_sub$SEX)
  rownames(annotation_sub) <- NULL
  SaveRDS(annotation_sub, paste0(OUT_DIR,"/selected_annotation_sub_",tt,".rds"))
  
  # entire GTEx
  mm <- model.matrix(data = annotation_sub, ~ SMTSISCH + SMGEBTCH + DTHHRDY + SEX)
  gtex_sub <- log2(gtex_sub+1); gc()
  gtex_sub <- apply(gtex_sub, 1, function(r) lm(r[as.numeric(rownames(mm))] ~ mm)$residuals); gc()
  SaveRDS(gtex_sub, paste0(OUT_DIR,"/log2_kown_batch_corrected_gtex_",tt,".rds"))
  

  #Reorganize annotation_sub
  gtex_sub <- t(gtex_sub)
  annotation_sub <- annotation_sub[annotation_sub$SAMPID %in% colnames(gtex_sub),]
  SaveRDS(annotation_sub, paste0(OUT_DIR,"/selected_annotation_sub_",tt,".rds"))
  
  gtex_final <- cbind(gtex_final,gtex_sub)
  annotation_final <- rbind(annotation_final,annotation_sub)
}

SaveRDS(gtex_final, paste0(OUT_DIR,"/log2_kown_batch_corrected_gtex.rds"))
SaveRDS(annotation_final, paste0(OUT_DIR,"/selected_annotation.rds"))







# Expression changes with age

DATA_DIR = "Outputs/3_GTExDataPrep/Subset_Data/Max_Subset"
OUTDIR = "Outputs/3_GTExDataPrep/Differential_Expression/Max_Subset"
INCL_70 = F

library(limma)
source("Scripts/functions.R")

metadata <- ReadRDS("Outputs/3_GTExDataPrep/metadata_restrictedaccess.rds")

files <- list.files(DATA_DIR, pattern = "sampled_data")
tissues <- sapply(files, function(x) strsplit(x, "_")[[1]][1])


for(tissue in unique(tissues)){
  
  # Read in expression for all age groups
  exprs <- sapply(paste0(DATA_DIR,"/",files[tissues == tissue]),
                  ReadRDS, simplify = F)
  exprs <- do.call(cbind, exprs)
  
  m <- subset(metadata, SAMPID %in% colnames(exprs))
  m$Age <- 25
  m$Age[m$AGE_GROUP == "30-39"] <- 35
  m$Age[m$AGE_GROUP == "40-49"] <- 45
  m$Age[m$AGE_GROUP == "50-59"] <- 55
  m$Age[m$AGE_GROUP == "60-69"] <- 65
  m$Age[m$AGE_GROUP == "70-79"] <- 70
  
  if(!INCL_70){
    m <- subset(m, AGE_GROUP != "70-79")
    exprs <- exprs[,m$SAMPID]
  }
  mm <- model.matrix(~ Age, data = m)
  
  fit <- limma::lmFit(exprs, mm)
  fit <- eBayes(fit)
  
  if(INCL_70){
    WriteRDS(fit, paste0(OUTDIR,"/",tissue,"_DE.rds"))
  } else{
    WriteRDS(fit, paste0(OUTDIR,"/",tissue,"_DE_no70.rds"))
  }
}

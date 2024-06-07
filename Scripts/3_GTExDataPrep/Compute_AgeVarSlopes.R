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
  vars <- lapply(exprs, function(m) apply(m, 1, var))
  vars <- do.call(cbind, vars)
  colnames(vars) <- sapply(files[tissues == tissue],
                           function(x) strsplit(x, "_")[[1]][2])
  
  if(INCL_70){
    lm_results <- t(apply(vars, 1, function(var){
      fit <- lm(data = data.frame("Var" = var,
                                  "Age" = c(25,35,45,55,65,70)[1:length(var)]),
                formula = Var ~ Age)
      coefs <- summary(fit)$coefficients
      return(c("Slope" = coefs["Age","Estimate"],
               "pval" = coefs["Age","Pr(>|t|)"]))
    }))
    lm_results <- cbind(lm_results,
                        "FDR" = p.adjust(lm_results[,"pval"], method = "fdr"))
    
    WriteRDS(lm_results, paste0(OUTDIR,"/",tissue,"_var_regression.rds"))
  } else{
    lm_results <- t(apply(vars, 1, function(var){
      fit <- lm(data = data.frame("Var" = var[1:5],
                                  "Age" = c(25,35,45,55,65)),
                formula = Var ~ Age)
      coefs <- summary(fit)$coefficients
      return(c("Slope" = coefs["Age","Estimate"],
               "pval" = coefs["Age","Pr(>|t|)"]))
    }))
    lm_results <- cbind(lm_results,
                        "FDR" = p.adjust(lm_results[,"pval"], method = "fdr"))
    
    WriteRDS(lm_results, paste0(OUTDIR,"/",tissue,"_var_regression_no70.rds"))
  }
}

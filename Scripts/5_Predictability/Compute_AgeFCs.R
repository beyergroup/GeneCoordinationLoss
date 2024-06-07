# Create linear model along age groups

args <- commandArgs(trailingOnly = T)
PATTERN = args[1]
COR_DIR = args[2]
OUTDIR = args[3]

source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")

if(COR_THRE == "well_predicted")
  cor_thre <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/correlation_thresholds.rds")


# Read in correlation files
files <- list.files(COR_DIR, pattern = PATTERN)
files <- files[grep("sampled_centered_correlations_spearman",files)]

tissues <- gsub("_maxsampled_centered_correlations_spearman.rds","",files)
tissues <- gsub("_old|_young","",tissues)

for(tissue in unique(tissues)){
  
  # Read in correlations for all age groups
  cors <- sapply(paste0(COR_DIR,"/",files[tissues == tissue]),
                 ReadRDS, simplify = T)
  colnames(cors) <- sapply(files[tissues == tissue],
                           function(x) strsplit(x, "_")[[1]][2])
  
  # Remove cases of NAs
  cors <- cors[rowSums(is.na(cors)) == 0,]
  
  # Remove poorly predicted genes
  if(COR_THRE == "well_predicted"){
    cors <- cors[rownames(cors) %in% well_predicted[[tissue]],]
  } else{
    cors <- cors[rowMeans(cors) > COR_THRE,]
  }
  
  # Linear model with all age groups
  lm_results <- t(apply(cors, 1, function(cor){
    fit <- lm(data = data.frame("Correlation" = cor,
                                "Age" = factor(c("old","young"),
                                               levels = c("young","old"))),
              formula = Correlation ~ Age)
    coefs <- summary(fit)$coefficients
    return(c("Slope" = coefs["Ageold","Estimate"],
             "pval" = coefs["Ageold","Pr(>|t|)"],
             "AdjustedSlope" = coefs["Ageold","Estimate"]/mean(cor)))
  }))
  lm_results <- cbind(lm_results,
                      "FDR" = p.adjust(lm_results[,"pval"], method = "fdr"))
  WriteRDS(lm_results, paste0(OUTDIR,"/",tissue,"_ageFC_",COR_THRE,".rds"))
  rm(lm_results,cors); gc()
}

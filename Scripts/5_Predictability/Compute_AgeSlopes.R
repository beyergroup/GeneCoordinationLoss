# Create linear model along age groups

args <- commandArgs(trailingOnly = T)
PATTERN = args[1]
COR_DIR = args[2]
INCL_70 = args[3]

source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")

COR_THRE = "well_predicted"

# Read in correlation files
files <- list.files(COR_DIR, pattern = PATTERN)
files <- files[grep("sampled_centered_correlations",files)]
# files <- files[-grep("young|old",files)]
files <- files[grep("spearman",files)]

# Read in well predicted genes
cor_thre <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/correlation_thresholds_spearman.rds")

tissues <- sapply(files, function(x) strsplit(x, "_")[[1]][1])

for(tissue in unique(tissues)){
  
  # Read in correlations for all age groups
  cors <- sapply(paste0(COR_DIR,"/",files[tissues == tissue]),
                 ReadRDS, simplify = T)
  colnames(cors) <- sapply(files[tissues == tissue],
                           function(x) strsplit(x, "_")[[1]][2])
  
  # Remove cases of NAs
  cors <- cors[rowSums(is.na(cors)) == 0,]
  
  if(!(tissue %in% names(cor_thre))){
    next
  }
  
  # Remove poorly predicted genes
  if(COR_THRE == "well_predicted"){
    cors <- cors[rowMeans(cors) > cor_thre[tissue],]
  } else{
    cors <- cors[rowMeans(cors) > COR_THRE,]
  }
  
  # Linear model with all age groups
  if(INCL_70){
    
    lm_results <- t(apply(cors, 1, function(cor){
      fit <- lm(data = data.frame("Correlation" = cor,
                                  "Age" = c(25,35,45,55,65,70)[1:length(cor)]),
                formula = Correlation ~ Age)
      coefs <- summary(fit)$coefficients
      return(c("Slope" = coefs["Age","Estimate"],
               "pval" = coefs["Age","Pr(>|t|)"],
               "AdjustedSlope" = coefs["Age","Estimate"]/mean(cor)))
      }))
    lm_results <- cbind(lm_results,
                        "FDR" = p.adjust(lm_results[,"pval"], method = "fdr"))
    WriteRDS(lm_results, paste0(COR_DIR,"/",tissue,
                                "_ageslope_",COR_THRE,".rds"))
    rm(lm_results,cors); gc()
    
  } else{
    
    lm_results <- t(apply(cors, 1, function(cor){
      fit <- lm(data = data.frame("Correlation" = cor[1:5],
                                  "Age" = c(25,35,45,55,65)),
                formula = Correlation ~ Age)
      coefs <- summary(fit)$coefficients
      return(c("Slope" = coefs["Age","Estimate"],
               "pval" = coefs["Age","Pr(>|t|)"],
               "AdjustedSlope" = coefs["Age","Estimate"]/mean(cor)))
    }))
    lm_results <- cbind(lm_results,
                        "FDR" = p.adjust(lm_results[,"pval"], method = "fdr"))
    WriteRDS(lm_results, paste0(COR_DIR,"/",tissue,
                                "_ageslope_",COR_THRE,"_no70.rds"))
    rm(lm_results,cors); gc()
  }
}

rm(files,tissues,tissue); gc()

# Define thresholds based on the background dist (for cross-tissue and per tissue)

args <- commandArgs(trailingOnly = T)
PATTERN = args[1]
REAL_DIR = args[2]
RAND_DIR = args[3]
MODE = args[4]

library(ggplot2)
source("Scripts/functions.R")

files <- list.files(REAL_DIR, pattern = PATTERN)
files <- files[grep("sampled_centered_correlations",files)]
if(MODE == "Spearman"){
  files <- files[grep("spearman",files)]
} else{
  files <- files[-grep("spearman",files)]
}

thresholds <- c()
predictable <- list()

for(file in files){
  
  tissue <- strsplit(file,"_")[[1]][1]
  message("\n",tissue)
  
  real <- ReadRDS(paste0(REAL_DIR,"/",file))
  rand <- ReadRDS(paste0(RAND_DIR,"/",file))
  rand <- rand[intersect(names(real),names(rand))]
  real <- real[names(rand)]
  
  # tail of 5%
  thre <- round(quantile(na.omit(rand), probs = 0.95),2)
  message("threshold = ",thre)
  thresholds <- c(thresholds,thre)
  
  message(sum(real > thre, na.rm = T), " correlations considered")
  
  predictable[[tissue]] <- names(which(real > thre))
}

names(thresholds) <- sapply(files, function(x) strsplit(x,"_")[[1]][1])

if(MODE == "Spearman"){
  WriteRDS(thresholds, paste0(REAL_DIR,"/correlation_thresholds_spearman.rds"))
  WriteRDS(predictable, paste0(REAL_DIR,"/well_predicted_genes_spearman.rds"))
} else{
  WriteRDS(thresholds, paste0(REAL_DIR,"/correlation_thresholds.rds"))
  WriteRDS(predictable, paste0(REAL_DIR,"/well_predicted_genes.rds"))
}

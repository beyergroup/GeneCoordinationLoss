# Define thresholds based on the background dist (for cross-tissue and per tissue)

args <- commandArgs(trailingOnly = T)
REAL_DIR = args[1]
RAND_DIR = args[2]

files <- list.files(REAL_DIR, pattern = "_sampled_centered_correlations.rds")

library(ggplot2)

thresholds <- c()
predictable <- list()

for(file in files){
  
  tissue <- strsplit(file,"_")[[1]][1]
  message("\n",tissue)
  
  real <- readRDS(paste0(REAL_DIR,file))
  rand <- readRDS(paste0(RAND_DIR,file))
  rand <- rand[intersect(names(real),names(rand))]
  real <- real[names(rand)]
  
  # tail of 5%
  thre <- round(quantile(na.omit(rand), probs = 0.95),2)
  message("threshold = ",thre)
  thresholds <- c(thresholds,thre)
  print(plot(density(na.omit(rand)), main = tissue))
  abline(v = thre)
  
  message(sum(real > thre, na.rm = T), " correlations considered")
  
  predictable[[tissue]] <- names(which(real > thre))
}

names(thresholds) <- sapply(files, function(x) strsplit(x,"_")[[1]][1])
saveRDS(thresholds, paste0(REAL_DIR,"/correlation_thresholds.rds"))
saveRDS(predictable, paste0(REAL_DIR,"/well_predicted_genes.rds"))

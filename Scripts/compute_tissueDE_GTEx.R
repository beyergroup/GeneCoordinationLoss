# DE with age in GTEx tissues

library(limma, lib.loc = "Resources/Rlibs/R-4.0.3/")

setwd("../")

files <- list.files("GTEx_Networks/Tissue_Networks/Outputs", pattern = "sampled_data.rds", full.names = T)
files <- files[-grep("CrossTissue",files)]

tissues <- unique(sapply(files, function(x) strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1]))

for(tissue in tissues){
  
  t <- readRDS(files[grep(tissue,files, fixed = T)])
  
  others <- sapply(files[-grep(tissue,files, fixed = T)], readRDS, simplify = F)
  others <- do.call(cbind,others)
  
  data <- cbind(t,others)
  metadata <- data.frame("Group" = c(rep(tissue,ncol(t)),rep("Others",ncol(others))))
  metadata$Group <- factor(metadata$Group, levels = c("Others",tissue))
  mm <- model.matrix(~ Group, data = metadata)
  
  fit <- limma::lmFit(data, mm)
  fit <- eBayes(fit)
  topTable(fit, coef = paste0("Group",tissue))
  
  saveRDS(fit, paste0("Outputs/GTEx/TissueDE/",tissue,"_sampled_tissueDE.rds"))
  rm(y,o,data,metadata,mm,fit); gc()
}


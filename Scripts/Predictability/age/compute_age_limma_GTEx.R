# DE with age in GTEx tissues

library(limma, lib.loc = "Resources/Rlibs/R-4.0.3/")

setwd("../")

# files <- list.files("GTEx_Networks/AgeTissue_Networks/Outputs", pattern = "sampled_data.rds", full.names = T)
files <- list.files("Outputs/Human_Network/stabsel/Predictability/AgeTissue", pattern = "sampled_delta.rds", full.names = T)

tissues <- unique(sapply(files, function(x) strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1]))

for(tissue in tissues){
  
  y <- cbind(readRDS(files[grep(paste0(tissue,"_20-29"),files, fixed = T)]),
             readRDS(files[grep(paste0(tissue,"_30-39"),files, fixed = T)]))
  
  o <- cbind(readRDS(files[grep(paste0(tissue,"_50-59"),files, fixed = T)]),
             readRDS(files[grep(paste0(tissue,"_60-69"),files, fixed = T)]))
  
  data <- cbind(y,o)
  metadata <- data.frame("AgeGroup" = c(rep("Young",ncol(y)),rep("Old",ncol(o))))
  metadata$AgeGroup <- factor(metadata$AgeGroup, levels = c("Young","Old"))
  
  mm <- model.matrix(~ AgeGroup, data = metadata)
  fit <- limma::lmFit(data, mm)
  fit <- eBayes(fit)
  print(topTable(fit, coef = "AgeGroupOld"))

  # coefficients <- sapply(rownames(data), function(gene){
  #   fit.data <- data.frame("Expression" = data[gene,], "AgeGroup" = metadata$AgeGroup)
  #   fit <- lm(Expression ~ AgeGroup, fit.data)
  #   return(summary(fit)$coefficients["AgeGroupOld",])
  # })
  # coefficients <- t(coefficients)
  # coefficients <- cbind(coefficients, "p.adj" = p.adjust(coefficients[,"Pr(>|t|)"], method = "fdr"))
  # 
  # message(tissue,"\n",min(coefficients[,"Pr(>|t|)"]))
  
  # saveRDS(coefficients, paste0("Outputs/GTEx/AgeDE/",tissue,"_sampled_ageLM.rds"))
  saveRDS(fit, paste0("Outputs/Human_Network/stabsel/Predictability/AgeTissue/",tissue,"_sampled_ageDP.rds"))
  # rm(y,o,data,metadata,coefficients); gc()
  
  # saveRDS(fit, paste0("Outputs/GTEx/AgeDE/",tissue,"_sampled_ageDE.rds"))
  rm(y,o,data,metadata,mm,fit); gc()
  
}

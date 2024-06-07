# Compute expression LFCs from one tissue to all others, for the whole genome

library(limma)
source("Scripts/5_Predictability/params.R")


# Load DESeq2 normalized data, un-centered
files <- paste0("Outputs/3_GTExDataPrep/Subset_Data/CrossAge/", ALL_TISSUES,
                "_sampled_data.rds")

for(tissue in ALL_TISSUES){
  
  t <- readRDS(files[grep(tissue, files, fixed = T)])
  
  others <- sapply(files[-grep(tissue, files, fixed = T)], readRDS,
                   simplify = F)
  others <- do.call(cbind,others)
  
  data <- cbind(t,others)
  metadata <- data.frame("Group" = c(rep(tissue, ncol(t)),
                                     rep("Others", ncol(others))))
  metadata$Group <- factor(metadata$Group, levels = c("Others", tissue))
  mm <- model.matrix(~ Group, data = metadata)
  
  fit <- limma::lmFit(data, mm)
  fit <- eBayes(fit)
  topTable(fit, coef = paste0("Group",tissue))
  
  saveRDS(fit,
          paste0("Outputs/5_Predictability/WellPredicted_TissueFilters/TissueDE/",
                 tissue, "_sampled_tissueDE.rds"))
  rm(t,others,data,metadata,mm,fit); gc()
}

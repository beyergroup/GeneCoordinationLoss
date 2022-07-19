# Compute correlation-based predictability LFCs

# FC = (1-corO)/(1-corY)
# LFC = log((1-corO)/(1-corY)) = log(1-corO) - log(1-corY)
# corrected LFC = 2*(log(1-corO) - log(1-corY))/(log(1-corO) + log(1-corY))

source("Scripts/functions.R")

well_predicted <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")

for(mode in c("cis","trans")){
  
  # files <- list.files("Outputs/5_Predictability",
  #                     pattern = paste0("_sampled_centered_",mode,
  #                                      "_correlations.rds"),
  #                     full.names = T)
  files <- list.files("Outputs/5_Predictability",
                      pattern = paste0("_sampled_centered_",mode,
                                       "_correlations_spearman.rds"),
                      full.names = T)
  old_files <- files[grep("old", files)]
  young_files <- files[grep("young", files)]

  for(i in 1:length(old_files)){
    
    tissue <- strsplit(tail(strsplit(old_files[i], split = "/")[[1]],1), "_")[[1]][2]
    
    old <- readRDS(old_files[i])
    young <- readRDS(young_files[i])
    
    genes <- intersect(names(old),names(young))
    genes <- intersect(genes,well_predicted[[tissue]])
    
    lfcs <- ComputeCorrelationLFC_adj(old[genes], young[genes])
    
    saveRDS(lfcs, paste0("Outputs/5_Predictability/",tissue,
                         "_",mode,"_corLFCs_YvsO_adj_spearman.rds"))
  }
}


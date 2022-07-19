# Compute NMSE-based predictability LFCs

# FC = NMSEo/NMSEy
# LFC = log(NMSEo/NMSEy) = log(NMSEo) - log(NMSEy)
# corrected LFC = 2*(log(NMSEo) - log(NMSEy))/(log(NMSEo) + log(NMSEy))

source("Scripts/functions.R")

well_predicted <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")

files <- list.files("Outputs/5_Predictability",
                    pattern = "_sampled_centered_NMSE.rds",
                    full.names = T)
old_files <- files[grep("old", files)]
young_files <- files[grep("young", files)]


for(i in 1:length(old_files)){
  
  tissue <- strsplit(tail(strsplit(old_files[i], split = "/")[[1]],1), "_")[[1]][1]
  
  old <- readRDS(old_files[i])
  young <- readRDS(young_files[i])
  
  genes <- intersect(names(old),names(young))
  genes <- intersect(genes,well_predicted[[tissue]])
  
  lfcs <- ComputeErrorLFC(old[genes], young[genes])
  
  saveRDS(lfcs, paste0("Outputs/5_Predictability/",tissue,"_errorLFCs_YvsO.rds"))
}

.libPaths("Resources/Rlibs/R-4.3.1")

files <- list.files("Outputs/5_Predictability/WellPredicted_TissueFilters",
                    pattern = "sampled_centered_correlations_spearman",
                    full.names = TRUE)

for(file in files){
  
  cor <- readRDS(file)
  cor <- sort(cor, decreasing = T)
  
  prefix <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  write.table(cor, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE,
              file = paste0("Outputs/5_Predictability/WellPredicted_TissueFilters/GSEA/",
                            prefix,".rnk"))
}

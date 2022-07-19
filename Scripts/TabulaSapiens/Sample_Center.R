# Reduce to same cell numbers

METHOD = "smartseq2"
NORM_DIR = "Outputs/Tabula_Sapiens/SCTransformNormalized/all_detected_genes/"
OUTDIR = "Outputs/Tabula_Sapiens/SampledSCTNorm/"

dir.create(OUTDIR)
norm_files <- list.files(NORM_DIR, pattern = METHOD)

# Get number of cells for each cell type
cell_numbers <- c()
for(file in norm_files){
  data <- readRDS(paste0(NORM_DIR,file))
  cell_numbers <- c(cell_numbers, ncol(data))
  rm(data); gc()
}
names(cell_numbers) <- norm_files
rm(norm_files); gc()


# Threshold for cell number
THRE = c("10X" = 1000, "smartseq2" = 200)[METHOD]
files <- names(which(cell_numbers >= THRE))


# Carry on with cell types passing threshold
for(file in files){
  
  data <- readRDS(paste0(NORM_DIR,file))
  
  # Sample to threshold cell numbers
  data <- data[,sample(1:ncol(data),THRE)]
  saveRDS(data, paste0(OUTDIR,gsub(".rds","_sampled_data.rds",file)))
  
  # Center per gene across cells
  means <- rowMeans(assay(data,"SCT"))
  mean_centered <- sweep(assay(data,"SCT"), 1, means, "-")
  saveRDS(mean_centered, paste0(OUTDIR,gsub(".rds","_sampled_meancentered_data.rds",file)))
  medians <- rowMedians(as.matrix(assay(data,"SCT")))
  median_centered <- sweep(assay(data,"SCT"), 1, medians, "-")
  saveRDS(median_centered, paste0(OUTDIR,gsub(".rds","_sampled_mediancentered_data.rds",file)))
  
  rm(data,means,mean_centered,medians,median_centered); gc()
}


# Create pseudobulk

INDIR = "../../Outputs/Tabula_Sapiens/SCTransformNormalized/all_detected_genes/"
OUTDIR = "../../Outputs/Tabula_Sapiens/Pseudobulk/"

dir.create(OUTDIR)

files <- list.files(INDIR)

pseudobulk <- list("Sum" = c(),
                   "Mean" = c(),
                   "Median" = c())
ncells <- c()

for(file in files){
  
  sce <- readRDS(paste0(INDIR,file))
  
  ncells <- c(ncells, ncol(sce))
  
  # sum
  pseudobulk[["Sum"]] <- cbind(pseudobulk[["Sum"]],
                               rowSums(assay(sce,"SCT")))
  # mean
  pseudobulk[["Mean"]] <- cbind(pseudobulk[["Mean"]],
                                rowMeans(assay(sce,"SCT")))
  # median
  pseudobulk[["Median"]] <- cbind(pseudobulk[["Median"]],
                                  rowMedians(as.matrix(assay(sce,"SCT"))))
}

rownames(pseudobulk$Median) <- rownames(pseudobulk$Mean)
pseudobulk <- lapply(pseudobulk, function(m){
  colnames(m) <- gsub(".rds","",files)
  return(m)
})

saveRDS(pseudobulk, paste0(OUTDIR,"pseudobulk_list.rds"))

# Re-normalize across cell types (CPM)

norm_pseudobulk <- lapply(pseudobulk, function(m){
  sf <- colSums(m)/1000000
  m <- sweep(m, 2, sf, "/")
  return(m)
})

saveRDS(norm_pseudobulk, paste0(OUTDIR,"normalized_pseudobulk_list.rds"))

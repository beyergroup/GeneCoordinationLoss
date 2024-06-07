
files <- list.files("GTEx_Networks/Tissue_Networks/Outputs",
                    pattern = "sampled_data.rds")


for(file in files){
  
  gtex <- readRDS(paste0("GTEx_Networks/Tissue_Networks/Outputs/",file))
  
  # split genes into 100 blocks
  gene.list <- list()
  for(i in 1:100){
    if((ceiling(nrow(gtex)/100)*i) <= nrow(gtex)){
      gene.list[[i]] <- rownames(gtex)[(1+(ceiling(nrow(gtex)/100)*(i-1))):(ceiling(nrow(gtex)/100)*i)]
    } else{
      gene.list[[i]] <- rownames(gtex)[(1+(ceiling(nrow(gtex)/100)*(i-1))):nrow(gtex)]
    }
  }
  
  # compute correlation by blocks
  cor.mat <- matrix(data = NA, nrow = nrow(gtex), ncol = nrow(gtex),
                    dimnames = list(rownames(gtex),rownames(gtex)))
  for(i in 1:100){
    for(j in i:100){
      current_cor.mat <- cor(t(gtex[gene.list[[i]],]), t(gtex[gene.list[[j]],]))
      cor.mat[rownames(current_cor.mat),colnames(current_cor.mat)] <- current_cor.mat
      rm(current_cor.mat); gc()
    }
  }
  
  saveRDS(cor.mat, paste0("Outputs/3_GTExDataPrep/",gsub("sampled_data", "correlation_mat",file)))
  rm(gtex, cor.mat, gene.list); gc()
}

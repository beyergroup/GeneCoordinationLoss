# Transform to quantile expression

WDIR = "3_TSDataPrep"
METHOD = "smartseq2"

.libPaths("Resources/Rlibs/R-4.0.3/")
library(SingleCellExperiment)
source("Scripts/functions.R")
source(paste0("Scripts/",WDIR,"/params.R"))

files <- list.files(paste0("Outputs/",WDIR,"/Normalized/Subset"),
                    pattern = paste0(METHOD,".rds"))

subset_cells <- list("10X" = ReadRDS("Outputs/3_TSDataPrep/subset_cells_10X.rds"),
                     "smartseq2" = ReadRDS("Outputs/3_TSDataPrep/subset_cells_smartseq2.rds"))
celltypes <- names(which(table(unlist(lapply(subset_cells, names))) == 2))
subset_cells <- lapply(subset_cells, function (l) l[names(l) %in% celltypes])

files <- files[grep(paste0(c("Merged_merged",celltypes), collapse = "|"),files)]

for(file in files){
  
  data <- ReadRDS(paste0("Outputs/",WDIR,"/Normalized/Subset/",file))
  data <- data[rowSums(data != 0) >= QCELL_THRE[METHOD],]
  message(nrow(data)," genes included")
  
  rnames <- c()
  quantile <- c()
  # compute quantiles
  for(i in 1:nrow(data)){
    tmp <- data[i,]
    if(sum(tmp != 0) >= QCELL_THRE[METHOD]){
      # compute quantiles excluding 0
      qs <- c(0, quantile(tmp[tmp != 0], probs = seq(0,1,.2))[1:5])
      if(length(unique(qs)) < length(qs)) # if 2 quantiles are the same
        next
      # convert into quantile
      tmp_q <- sapply(tmp, function(x) sum(x >= qs)-1)
      quantile <- rbind(quantile, tmp_q)
      rnames <- c(rnames,rownames(data)[i])
    }
  }
  rm(tmp,qs,tmp_q); gc()
  
  rownames(quantile) <- rnames
  saveRDS(quantile, paste0("Outputs/",WDIR,"/Normalized/Subset/",
                           gsub(".rds","_quantile.rds",file)))
  rm(data,quantile,rnames); gc()
}

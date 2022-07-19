# .libPaths("../../Resources/Rlibs/R-4.0.3/previous_Seurat")
# library(Seurat)
.libPaths("../../Resources/Rlibs/R-4.0.3/Seurat-3-2-3")
library(Seurat)

INDIR = "../../Data/Tabula_Sapiens/tabula_sapiens"
OUTDIR = "../../Outputs/Tabula_Sapiens"

# read the full data from .rds files
files <- list.files(INDIR, full.names = T)

# create temporary directory where data is separated by compartment and tissue
OUTDIR <- paste0(OUTDIR,"/tmp")
dir.create(OUTDIR)

for(file in files){
  
  seu <- readRDS(file)
  
  # adjust colnames of nCounts and nFeatures
  if(length(grep("nCount", colnames(seu@meta.data))) == 1){
    colnames(seu@meta.data)[grep("nCount", colnames(seu@meta.data))] <- "nCount"
  } else{
    stop("Not sure which column to take for nCount")
  }
  if(length(grep("nFeature", colnames(seu@meta.data))) == 1){
    colnames(seu@meta.data)[grep("nFeature", colnames(seu@meta.data))] <- "nFeature"
  } else{
    stop("Not sure which column to take for nFeature")
  }
  
  seu$Age <- as.numeric(sapply(as.character(seu$development_stage),
                               function(x) strsplit(x,"-")[[1]][1]))
  Idents(seu) <- "tissue"
  
  tissues <- as.character(unique(seu$tissue))
  for(t in tissues){
    
    # manually remove key bullshit bc it crashes the merging
    seu@assays$RNA@key <- t
    
    saveRDS(subset(seu, idents = t),
            paste0(OUTDIR,"/",t,"_",as.character(unique(seu$compartment)),".rds"))
  }
  rm(seu,tissues,t); gc()
}
rm(file,files); gc()


# read from temporary directory and join compartments of the same tissue

# detach("package:Seurat", unload=TRUE)
# detach("package:SeuratObject", unload=TRUE)
# .libPaths("../../Resources/Rlibs/R-4.0.3/")
# library(Seurat)

files <- list.files(OUTDIR)
tissues <- unique(sapply(files, function(x) strsplit(x, "_")[[1]][1]))

for(t in tissues){
  
  seu_list <- sapply(paste0(OUTDIR,"/",files[grep(t,files)]), readRDS)
  
  # # arrange metadata annotation
  # seu_list <- lapply(seu_list, function(s){
  #   
  #   # manually remove key bullshit bc it crashes the merging
  #   s@assays$RNA@key <- t
  #   
  #   # adjust colnames of nCounts and nFeatures
  #   if(length(grep("nCount", colnames(s@meta.data))) == 1){
  #     colnames(s@meta.data)[grep("nCount", colnames(s@meta.data))] <- "nCount"
  #   } else{
  #     stop("Not sure which column to take for nCount")
  #   }
  #   if(length(grep("nFeature", colnames(s@meta.data))) == 1){
  #     colnames(s@meta.data)[grep("nFeature", colnames(s@meta.data))] <- "nFeature"
  #   } else{
  #     stop("Not sure which column to take for nFeature")
  #   }
  #   
  #   return(s)
  # })
  
  if(length(seu_list) > 2){
    seu <- merge(x = seu_list[[1]], y = seu_list[-1], merge.data = TRUE)
  } else{
    seu <- seu_list[[1]]
  }
  rm(seu_list); gc()
  
  saveRDS(seu,
          paste0(gsub(OUTDIR, pattern = "tmp",
                      replacement = "pertissue"),
                 "/",t,".rds"))
  rm(seu); gc()
}

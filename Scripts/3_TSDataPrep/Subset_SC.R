# Subset single cell data to right cell types and numbers

WDIR = "3_TSDataPrep"
METHOD = "10X"

source("Scripts/functions.R")
source(paste0("Scripts/",WDIR,"/params.R"))

files <- list.files(paste0("Outputs/",WDIR,"/Normalized"),
                    pattern = paste0(METHOD,".rds"))


for(file in files){
  
  # read in
  sce <- ReadRDS(paste0("Outputs/",WDIR,"/Normalized/",file))
  
  # check if enough cells of cell type have a rich transcriptome
  if(sum(sce$n_genes > NDETGENES) < NCELL_THRE_SC[METHOD]){
    next
  }
  
  # if yes, restrict to rich transcriptome cells
  sce <- sce[,sce$n_genes > NDETGENES]

  # subset to same cell numbers for all cell types
  data <- assay(sce,"SCT")[,sample(1:ncol(sce),NCELL_THRE_SC[METHOD])]

  # save
  WriteRDS(data, paste0("Outputs/3_TSDataPrep/Normalized/Subset/",file))

  rm(sce,data); gc()
}

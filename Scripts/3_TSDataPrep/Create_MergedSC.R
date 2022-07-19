# Create cross-tissue

WDIR = "3_TSDataPrep"

source("Scripts/functions.R")
source("Scripts/3_TSDataPrep/params.R")
library(SingleCellExperiment)
library(Seurat)
library(purrr)

# check cell types that passed filters from 10X and Smartseq2
subset_cells <- list("10X" = ReadRDS("Outputs/3_TSDataPrep/subset_cells_10X.rds"),
                     "smartseq2" = ReadRDS("Outputs/3_TSDataPrep/subset_cells_smartseq2.rds"))
celltypes <- names(which(table(unlist(lapply(subset_cells, names))) == 2))
subset_cells <- lapply(subset_cells, function (l) l[names(l) %in% celltypes])

# pick up files prior to normalization
files <- list.files(paste0("Outputs/",WDIR,"/Filters/"),
                    pattern = paste0(celltypes, collapse = "|"))

# merge per sequencing method
for(method in c("smartseq2","10X")){
  sce_list <- sapply(files[grep(method,files)],
                     function(f) ReadRDS(paste0("Outputs/",WDIR,
                                                "/Filters/",f)))
  # convert to Seurat
  seu_list <- lapply(sce_list, function(sce){
    CreateSeuratObject(counts = as.matrix(assay(sce,"decontXcounts")),
                       meta.data = as.data.frame(sce@colData))
  })
  rm(sce_list); gc()
  seu <- reduce(seu_list, merge)
  rm(seu_list); gc()
  WriteRDS(seu, paste0("Outputs/",WDIR,
                       "/Filters/alldetected_genes_Merged_merged_",
                       method,"_seu.rds"))
  
  # subset per cell type
  cells <- c()
  ncells <- floor(NCELL_THRE_SC[[method]]/
                    length(names(subset_cells[[method]])))
  for(ct in names(subset_cells[[method]])){
    cells <- c(cells, sample(colnames(seu)[paste(seu$organ_tissue,
                                                 gsub(" ","_",seu$cell_ontology_class),
                                                 sep = "_") == ct], ncells))
  }
  seu <- seu[,cells]
  
  # normalize
  seu <- SCTransform(seu, method = "glmGamPoi", vst.flavor = "v2")
  WriteRDS(seu@assays$SCT@data, paste0("Outputs/",WDIR,
                        "/Normalized/Subset/alldetected_genes_Merged_merged_",
                        method,".rds"))
}

# get hepatocytes from 10X
sce_liver <- readRDS("../../Outputs/Tabula_Sapiens/pertissue/10X/Liver.rds")
sce_liver <- sce_liver[, sce_liver@colData$cell_ontology_class == "hepatocyte", drop=TRUE]

# get colon epithelial cells from 10X
sce_colon <- readRDS("../../Outputs/Tabula_Sapiens/pertissue/10X/Large_Intestine.rds")
sce_colon <- sce_colon[, sce_colon@colData$compartment == "endothelial", drop = TRUE]

saveRDS(sce_liver,"/cellnet/Luise/Carolina/TabulaSapiens_liver.rds")
saveRDS(sce_colon,"/cellnet/Luise/Carolina/TabulaSapiens_colon.rds")

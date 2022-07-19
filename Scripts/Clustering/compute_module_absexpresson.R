# Compute module cross-tissue absolute expression

source("Scripts/functions.R")

# # take cross tissue samples
# gtex <- readRDS("GTEx_Networks/Tissue_Networks/Outputs/CrossTissue_sampled_data.rds")
# samples <- colnames(gtex)
# rm(gtex); gc()
# # read in lib-size normalized expr. values
# gtex_norm <- readRDS("Data/GTEx/DESeq2_normalized_gtex.rds")
# gtex_norm <- gtex_norm[,samples]
# rm(samples); gc()
# # convert to gene symbol
# gtex_norm <- DataENSGToSymbol(gtex_norm, remove_dup = T)
# gtex_norm <- log2(1+gtex_norm)
# 
# for(NET in c("stabsel","stabsel_pcclasso")){
#   
#   for(METHOD in c("greedy","rwalk")){
#     
#     # module membership
#     membership <- readRDS(paste0("Outputs/Human_Network/",NET,
#                                  "/Topology/Modules/membership_absweights_",
#                                  METHOD,"_iterreclustering.rds"))
#     modules <- names(table(membership)[table(membership) >= 10])
#     
#     abs_expression <- c()
#     
#     for(module in modules){
#       genes <- names(which(membership == module))
#       abs_expression <- c(abs_expression,
#         mean(gtex_norm[intersect(genes,rownames(gtex_norm)),]))
#     }
#     names(abs_expression) <- modules
#     
#     saveRDS(abs_expression, paste0("Outputs/Human_Network/",NET,
#       "/Topology/Modules/TissueDE/CrossTissue_mean_absexpression_",METHOD,".rds"))
#   }
# }


# Compute tissue-specific absolute expression

tissue_files <- list.files("GTEx_Networks/Tissue_Networks/Outputs",
                           pattern = "sampled_data.rds", full.names = T)

gtex_norm <- readRDS("Data/GTEx/DESeq2_normalized_gtex.rds")

for(file in tissue_files){
  
  # read in tissue-specific GTEx subset for sample identification
  gtex <- readRDS(file)
  samples <- colnames(gtex)
  rm(gtex); gc()
  
  # subset big GTEx
  gtex <- gtex_norm[,samples]
  rm(samples); gc()
  # convert to gene symbol
  gtex <- DataENSGToSymbol(gtex, remove_dup = T)
  gtex <- log2(1+gtex)
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  for(NET in c("stabsel","stabsel_pcclasso")){
    
    for(METHOD in c("greedy","rwalk")){
      
      # module membership
      membership <- readRDS(paste0("Outputs/Human_Network/",NET,
                                   "/Topology/Modules/membership_absweights_",
                                   METHOD,"_iterreclustering.rds"))
      modules <- names(table(membership)[table(membership) >= 10])
      
      abs_expression <- c()
      
      for(module in modules){
        genes <- names(which(membership == module))
        abs_expression <- c(abs_expression,
                            mean(gtex[intersect(genes,rownames(gtex)),]))
      }
      names(abs_expression) <- modules
      
      saveRDS(abs_expression, paste0("Outputs/Human_Network/",NET,
        "/Topology/Modules/TissueDE/",tissue,"_mean_absexpression_",METHOD,".rds"))
    }
  }
}

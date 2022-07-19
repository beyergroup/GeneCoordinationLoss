# Run GO enrichment analysis

args <- commandArgs(trailingOnly = T)
NET = args[1]
WDIR = args[2]
PATTERN = args[3]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")
source(paste0("Scripts/",WDIR,"/params.R"))
library(ggplot2)
library(ggpubr)

if(PATTERN == "all"){
  files <- list.files(paste0("Outputs/",WDIR))
} else{
  # add error for when PATTERN is not among files
  files <- list.files(paste0("Outputs/",WDIR),
                      pattern = PATTERN)
  files <- files[grep(NET,files)]
}
# files <- files[grep("_membership.rds",files)]
files <- files[grep("_membership_reclustering_redecomposing_250_adjusted.rds",files)]

for(file in files){
  
  # read in module membership
  membership <- ReadRDS(paste0("Outputs/",WDIR,"/",file))
  
  # limit to big enough modules
  modules <- names(table(membership)[table(membership) >= 10])
  
  GObp.list <- list()
  GOmf.list <- list()
  GOcc.list <- list()
  
  for(module in modules){
    
    genes <- names(which(membership == module))
    
    GObp.list[[module]] <- GetGOEnrich(genes, names(membership),
                                       "BP", enrich_cutoff = 1,
                                       algorithm = "weight01")
    GOmf.list[[module]] <- GetGOEnrich(genes, names(membership),
                                       "MF", enrich_cutoff = 1,
                                       algorithm = "weight01")
    GOcc.list[[module]] <- GetGOEnrich(genes, names(membership),
                                       "CC", enrich_cutoff = 1,
                                       algorithm = "weight01")
  }
  save(GObp.list,GOmf.list,GOcc.list,
       file = paste0("Outputs/",WDIR,"/",
                     gsub(".rds","_GO_weight01.RData",file)))
}

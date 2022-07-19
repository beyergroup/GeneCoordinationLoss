# Define clusters (cut tree at chosen height)

args <- commandArgs(trailingOnly = TRUE)
NET = args[1]
WDIR = args[2]
PATTERN = args[3]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")

if(PATTERN == "all"){
  files <- list.files(paste0("Outputs/",WDIR))
} else{
  # add error for when PATTERN is not among files
  files <- list.files(paste0("Outputs/",WDIR),
                      pattern = PATTERN)
  files <- files[grep(NET,files)]
}
files <- files[grep("_height_eigen_decision.rds",files)]


for(file in files){
  
  # load decisions
  decision <- ReadRDS(paste0("Outputs/",WDIR,"/",file))
  
  # load clustering with chosen number of eigenvectors
  hc <- ReadRDS(paste0("Outputs/",WDIR,"/",
                       gsub("height_eigen_decision",
                            paste0("evectors_evalues_complete_link_clustering_",
                                   as.character(decision$Eigenvectors),
                                   "_evectors"),file)))
  
  # cut at chosen height
  membership <- cutree(tree = hc, h = decision$Height)
  
  # save
  WriteRDS(membership, paste0("Outputs/",WDIR,"/",
                              gsub("height_eigen_decision","membership",file)))
}

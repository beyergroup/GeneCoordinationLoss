# Grab smoothed LFCs and average over modules

.libPaths("Resources/Rlibs/R-4.0.3/")
library(devtools)
library(igraph)

LFC_DIR = "/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel/Modules/Module_Activity/SmoothedLFCs/"
OUT_DIR = "/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel/Modules/Module_Activity/GlobalSmoothedActivity/"
MEMBERSHIP_FILE = "Outputs/Human_Network/stabsel/Modules/old/Adjacency_weightnone_membership.rds"

TS = T

dir.create(OUT_DIR)
if(TS){
  dir.create(OUT_DIR,"Tabula_Sapiens")
}

files <- list.files(LFC_DIR)

# Get module membership -------------------------------------------------------

membership <- readRDS(MEMBERSHIP_FILE)
modules <- names(table(membership)[table(membership) >= 10])

# Get smoothed values ---------------------------------------------------------

for(file in files){
  
  smoothed <- readRDS(paste0(LFC_DIR,file))
  
  # compute module activity
  a <- matrix(nrow = length(modules), ncol = ncol(smoothed),
              dimnames = list("Modules" = modules, "Tissues" = colnames(smoothed)))
  for(module in modules){
    genes <- intersect(names(which(membership == module)),rownames(smoothed))
    a[module,] <- apply(smoothed, 2, function(x) mean(x[genes]))
  }
  saveRDS(a, paste0(OUT_DIR,gsub(".rds","_relative_activity_GTEx.rds",file)))
}



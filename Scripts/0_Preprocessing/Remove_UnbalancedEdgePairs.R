# Remove smaller edge of unbalanced edge pairs that should be uni-directional

args <- commandArgs(trailingOnly = T)
NET = args[1]
WDIR = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")
source(paste0("Scripts/",WDIR,"/params.R"))

message("Edge weight / sum of pair threshold = ",EDGE_WEIGHTS_RATIO_THRE)

# Load network ----------------------------------------------------------------
net <- ReadRDS(paste0("Data/Networks/",NET,"_network_Hs.rds"))
net <- as.matrix(net)
net <- MatrixToSquare(net)

# Remove residual edges -------------------------------------------------------
net <- RemoveResEdges(net, threshold = EDGE_WEIGHTS_RATIO_THRE)

# Save output -----------------------------------------------------------------
WriteRDS(net, paste0("Outputs/",WDIR,"/",NET,"_filtered_network_Hs.rds"))

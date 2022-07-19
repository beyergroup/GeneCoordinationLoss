# Get Human networks

args <- commandArgs(trailingOnly = T)
NETLIBS = args[1]
NET = args[2]
WDIR = args[3]

source("Scripts/functions.R")

net <- ReadNetwork(STABSEL, NETLIBS)
saveRDS(net, paste0("Outputs/",WDIR,"/stabsel_network_Hs.rds"))
rm(net); gc()


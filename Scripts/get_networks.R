# Get Human networks

source("/data/public/xwu2/Rscripts/network_funclib.r")

## Stability selection network
NET_FILE = "/data/public/xwu2/networks/infer_CCLERSEQ_01/Network_Whole_CCLERSEQ_01.Rout"
network <- loadNetworkWithFilteringForSignificantPredictors(NET_FILE, pValCutoff=-1, localGeneCutoff = 50)
net <- network$G
rownames(net) <- network$responseGenes
colnames(net) <- network$predictors
rm(network); gc()
net <- net[,-1]
net <- net[rowSums(net) != 0,]
net <- net[,colSums(net) != 0]
saveRDS(net, "Outputs/Human_Network/stabsel_network_Hs.rds")
rm(NET_FILE,net); gc()

## No partial correlations network
NET_FILE = "/data/public/xwu2/networks/infer_CCLERSEQ_01/Network_Whole_CCLERSEQ_01.pcc_lasso_regress.regress.Rout"
network <- loadNetworkWithFilteringForSignificantPredictors(NET_FILE, pValCutoff=-1, localGeneCutoff = 50)
net <- network$G
rownames(net) <- network$responseGenes
colnames(net) <- network$predictors
rm(network); gc()
net <- net[,-1]
net <- net[rowSums(net) != 0,]
net <- net[,colSums(net) != 0]
saveRDS(net, "Outputs/Human_Network/stabsel_pcclasso_network_Hs.rds")

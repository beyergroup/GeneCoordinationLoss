# Compare in-going edges to out-going edges

net_pcc <- readRDS("/data/public/adesous1/GeneCorrelation/Data/Networks/Human/stabsel_pcclasso_network_Hs.rds")
net_pcc <- as.matrix(net_pcc)

e <- apply(which(upper.tri(net_pcc), arr.ind = T), 1, function(x) c(net_pcc[x[1],x[2]], net_pcc[x[2],x[1]]))
e <- t(e)
e <- e[rowSums(e) != 0,]

png("Plots/Human_Network/stabsel_pcclasso/Topology/in_out_edges_scatter.png")
# plot(x = e[,1], y = e[,2], xlab = "", ylab = "", main = "In- vs out-going edges")
smoothScatter(x = e[,1], y = e[,2], xlab = "", ylab = "", main = "In- vs out-going edges")
dev.off()

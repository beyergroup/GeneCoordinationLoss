# Compare in-going edges to out-going edges

WDIR = "1_Topology/"
args <- commandArgs(trailingOnly = T)
NET = args[1]

net_pcc <- readRDS("/data/public/adesous1/GeneCorrelation/Data/Networks/Human/stabsel_pcclasso_network_Hs.rds")
net_pcc <- as.matrix(net_pcc)

e <- apply(which(upper.tri(net_pcc), arr.ind = T), 1, function(x) c(net_pcc[x[1],x[2]], net_pcc[x[2],x[1]]))
e <- t(e)
e <- e[rowSums(e) != 0,]

png("Plots/Human_Network/stabsel_pcclasso/Topology/in_out_edges_scatter.png")
# plot(x = e[,1], y = e[,2], xlab = "", ylab = "", main = "In- vs out-going edges")
smoothScatter(x = e[,1], y = e[,2], xlab = "", ylab = "", main = "In- vs out-going edges")
dev.off()

ratios <- apply(e, 1,
  function(x) min(abs(x))/sum(abs(x)))

plot.data <- data.frame("Edge1" = e[,1], "Edge2" = e[,2], "Ratios" = ratios)

library(ggplot2, lib.loc = "/data/public/adesous1/GeneCorrelation/Resources/Rlibs/R-4.0.3/")
ggplot(subset(plot.data, Ratios > 0.1)) + geom_point(aes(x = Edge1, y = Edge2, color = Ratios))

n <- readRDS("/data/public/adesous1/GeneCorrelation/Data/Networks/Human/stabsel_network_Hs.rds")
n <- as.matrix(n)
net <- matrix(data = 0,
              nrow = length(union(colnames(n),rownames(n))),
              ncol = length(union(colnames(n),rownames(n))),
              dimnames = list(union(colnames(net),rownames(net)),
                              union(colnames(net),rownames(net))))
net[rownames(n),colnames(n)] <- n
rm(n): gc()

e <- apply(which(upper.tri(net), arr.ind = T), 1, function(x) c(net[x[1],x[2]], net[x[2],x[1]]))
e <- t(e)
e <- e[rowSums(e) != 0,]
table(rowSums(e == 0))
e <- e[rowSums(e == 0) == 0,]

png("Plots/Human_Network/stabsel/Topology/in_out_edges_nozero.png")
# plot(x = e[,1], y = e[,2], xlab = "", ylab = "", main = "In- vs out-going edges")
smoothScatter(x = e[,1], y = e[,2], xlab = "", ylab = "", main = "In- vs out-going edges")
dev.off()





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




# confirm edges in pcc network are a subset of the edges in the full network
net <- net[rownames(net_pcc),colnames(net_pcc)]
edges_pcc <- which(net_pcc != 0, arr.ind = T)
table(apply(edges_pcc, 1, function(x) net[x[1],x[2]] != 0)) # all edges of the pcc network exist in the full network
table(apply(edges_pcc, 1, function(x) net[x[1],x[2]] == net_pcc[x[1],x[2]])) # 3526 edges have the same weight
# what about the remaining 57112?
diff_edges <- which(apply(edges_pcc, 1, function(x) net[x[1],x[2]] != net_pcc[x[1],x[2]]))
diff_weights <- t(apply(edges_pcc[diff_edges,], 1, function(x) c("full" = net[x[1],x[2]],
                                                                 "pcc" = net_pcc[x[1],x[2]])))
plot(x = diff_weights[,"full"], y = diff_weights[,"pcc"])
abline(a = 0, b = 1)
# edges that are different get larger in pcc network (makes sense to account for less predictors)
# there are some exceptions, which probably correspond to the edges that become 'negligible'
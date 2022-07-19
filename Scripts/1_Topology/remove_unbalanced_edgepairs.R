source("Scripts/functions.R")

# network without indirect effects

net_pcc <- readRDS("/data/public/adesous1/GeneCorrelation/Data/Networks/Human/stabsel_pcclasso_network_Hs.rds")
net_pcc <- as.matrix(net_pcc)

ind <- which(upper.tri(net_pcc), arr.ind = T)

ratios <- apply(ind, 1,
  function(x) min(abs(c(net_pcc[x[1],x[2]],net_pcc[x[2],x[1]])))/sum(abs(c(net_pcc[x[1],x[2]],net_pcc[x[2],x[1]]))))

THRE = 0.1 # ratios below this threshold get turned into directed edges

to_remove <- ind[(ratios < THRE) & (!is.na(ratios)),] # 1244 out of 30319 edges removed
to_remove <- t(apply(to_remove, 1, function(x)
  switch(which.min(c(abs(net_pcc[x[1],x[2]]),abs(net_pcc[x[2],x[1]]))),
         c(x[1],x[2]),
         c(x[2],x[1]))))

net_pcc[to_remove] <- 0
saveRDS(net_pcc, "/data/public/adesous1/GeneCorrelation/Data/Networks/Human/stabsel_pcclasso_network_Hs_filtered.rds")

rm(ind,to_remove,ratios); gc()

e <- apply(which(upper.tri(net_pcc), arr.ind = T), 1, function(x) c(net_pcc[x[1],x[2]], net_pcc[x[2],x[1]]))
e <- t(e)
e <- e[rowSums(e == 0) == 0,]

smoothScatter(x = e[,1], y = e[,2], xlab = "", ylab = "", main = "In- vs out-going edges")
plot(x = e[,1], y = e[,2], xlab = "", ylab = "", main = "In- vs out-going edges")



## full network

n <- readRDS("/data/public/adesous1/GeneCorrelation/Data/Networks/Human/stabsel_network_Hs.rds")
n <- as.matrix(n)
net <- matrix(data = 0,
              nrow = length(union(colnames(n),rownames(n))),
              ncol = length(union(colnames(n),rownames(n))),
              dimnames = list(union(colnames(n),rownames(n)),
                              union(colnames(n),rownames(n))))
net[rownames(n),colnames(n)] <- n
rm(n); gc()

net <- RemoveResEdges(net, threshold = 0.1)
saveRDS(net, "/data/public/adesous1/GeneCorrelation/Data/Networks/Human/stabsel_network_Hs_filtered.rds")

e <- apply(which(upper.tri(net), arr.ind = T), 1, function(x) c(net[x[1],x[2]], net[x[2],x[1]]))
e <- t(e)
e <- e[rowSums(e == 0) == 0,]

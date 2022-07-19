# Compare in-going edges to out-going edges

args <- commandArgs(trailingOnly = T)
NET = args[1]
WDIR = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(ggplot2)
source("Scripts/functions.R")

# Load network ----------------------------------------------------------------
net <- ReadRDS(paste0("Data/Networks/",NET,"_network_Hs.rds"))
net <- as.matrix(net)


# Collect values for edge pairs (node i -> node j and opposite direction) -----
e <- apply(which(upper.tri(net), arr.ind = T), 1,
           function(x) c(net[x[1],x[2]], net[x[2],x[1]]))
e <- t(e)
e <- e[rowSums(e) != 0,]

# Plot
png(paste0("Plots/",WDIR,"/",NET,"_in_out_edges_scatter.png"))
smoothScatter(x = e[,1], y = e[,2], xlab = "", ylab = "",
              main = "In- vs out-going edges")
dev.off()


# Compute ratio between smallest edge and sum of both -------------------------
ratios <- apply(e, 1, function(x) min(abs(x))/sum(abs(x)))

# Plot
plot.data <- data.frame("Edge1" = e[,1],
                        "Edge2" = e[,2],
                        "Ratios" = ratios)

png(paste0("Plots/",WDIR,"/",NET,"_in_out_edges_ratios.png"))
ggplot(subset(plot.data, Ratios > 0.1)) +
  geom_point(aes(x = Edge1, y = Edge2, color = Ratios))
dev.off()


# Restrict to pairs where none of the weights is 0 ----------------------------
e <- e[rowSums(e == 0) == 0,]

# Plot
png(paste0("Plots/",WDIR,"/",NET,"_in_out_edges_nozero.png"))
smoothScatter(x = e[,1], y = e[,2], xlab = "", ylab = "",
              main = "In- vs out-going edges")
dev.off()





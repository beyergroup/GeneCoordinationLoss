membership <- readRDS("/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel/Modules/old/Adjacency_weightnone_membership.rds")

sizes <- as.numeric(table(membership))
names(sizes) <- names(table(membership))
sizes <- sizes[sizes >= 10]
sizes <- sort(sizes, decreasing = T)

pdf("Plots/Human_Network/stabsel/Modules/module_sizes.pdf", height = 5, width = 5)
barplot(sizes, xlab = "Module ID", ylab = "Module size")
dev.off()

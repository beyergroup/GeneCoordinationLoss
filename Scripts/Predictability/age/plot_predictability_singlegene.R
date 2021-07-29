# Plot delta across age groups for individual genes

files <- list.files("Outputs/Human_Network/stabsel/Predictability/AgeTissue",
                    pattern = "_sampled_delta.rds", full.names = T)

deltas <- sapply(files, readRDS, simplify = F)


gene = "LINC01927"

plot.data <- lapply(deltas, function(m) m[gene,])
plot.data <- melt(plot.data)
plot.data$Tissue <- sapply(plot.data$L1, function(x) strsplit(tail(strsplit(x,"/")[[1]],1), "_")[[1]][1])
plot.data$Age <- sapply(plot.data$L1, function(x) strsplit(tail(strsplit(x,"/")[[1]],1), "_")[[1]][2])

pdf(paste0("Plots/Human_Network/stabsel/Predictability/delta_DE_age/",gene,"_deltas.pdf"),
    height = 9, width = 15)
ggplot(plot.data) +
  geom_violin(aes(x = Age, y = value, fill = Age), size = 1, draw_quantiles = 0.5) +
  scale_fill_viridis_d(option = "magma", begin = .3, end = .7) +
  facet_wrap(~ Tissue) + ggtitle(gene) +
  ylab("Deviation") + xlab("Age") + guides(fill = F) +
  theme(text = element_text(size = 18), title = element_text(size = 18))
dev.off()

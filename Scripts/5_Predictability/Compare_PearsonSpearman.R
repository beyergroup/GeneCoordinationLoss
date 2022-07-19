# Compare LFCs based on Pearson and Spearman correlations

pearson_LFCs <- list.files("Outputs/5_Predictability",
                           pattern = paste0("corLFCs_YvsO_adj.rds"),
                           full.names = T)
spearman_LFCs <- list.files("Outputs/5_Predictability",
                            pattern = paste0("corLFCs_YvsO_adj_spearman.rds"),
                            full.names = T)

tissues <- unique(sapply(pearson_LFCs, function(x)
  strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1]))

plot.data <- data.frame()

for(tissue in tissues){
  
  pearson <- readRDS(pearson_LFCs[grep(tissue,pearson_LFCs,fixed = T)])
  spearman <- readRDS(spearman_LFCs[grep(tissue,spearman_LFCs,fixed = T)])
  
  plot.data <- rbind.data.frame(plot.data,
                                data.frame("Pearson" = pearson,
                                           "Spearman" = spearman,
                                           "Tissue" = tissue))
}

pdf("Plots/5_Predictability/spearman_pearson_LFCs.pdf",
    height = 10, width = 16)
ggplot(plot.data) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(aes(x = Pearson, y = Spearman),
             alpha = .3, size = .5) +
  facet_wrap(~ Tissue) +
  xlab("Pearson-based LFC") + ylab("Spearman-based LFC") +
  theme(text = element_text(size = 20))
dev.off()


# values pretty comparable

dif.data <- subset(plot.data, abs(plot.data$Pearson - plot.data$Spearman) > .5)
# previous hits that are lost


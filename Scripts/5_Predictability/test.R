
source("Scripts/functions.R")
library(ggplot2)

well_predicted_genes <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")

# old LFCs
old_LFCs <- list.files("Outputs/5_Predictability",
                       pattern = "_corLFCs_YvsO.rds",
                       full.names = T)
old_LFCs <- old_LFCs[-grep("cis|trans",old_LFCs)]

# adjusted LFCs
new_LFCs <- list.files("Outputs/5_Predictability",
                       pattern = "_corLFCs_YvsO_adj.rds",
                       full.names = T)

files <- list.files("Outputs/5_Predictability",
                    pattern = "_sampled_centered_correlations.rds",
                    full.names = T)
old_files <- files[grep("old", files)]
young_files <- files[grep("young", files)]
rm(files); gc()

plot.data <- data.frame()

for(i in 1:length(old_files)){
  
  tissue <- strsplit(tail(strsplit(old_files[i], split = "/")[[1]],1), "_")[[1]][2]
  
  old <- readRDS(old_files[i])
  young <- readRDS(young_files[i])
  
  oldLFC <- readRDS(old_LFCs[i])
  newLFC <- readRDS(new_LFCs[i])
  
  common <- intersect(intersect(names(young),names(old)),
                      intersect(names(oldLFC),names(newLFC)))
  common <- intersect(intersect(names(young),names(old)),
                      well_predicted_genes[[tissue]])

  message(tissue)
  print(table(old[common] < young[common],
              oldLFC[common] > 0))

  problem_genes <- names(which((old[common] < young[common]) &
                                 (oldLFC[common] < 0)))
  
  plot.data <- rbind.data.frame(plot.data,
                                data.frame("OldCor" = old[common],
                                           "YoungCor" = young[common],
                                           "PrevLFCs" = oldLFC[common],
                                           "AdjLFCs" = newLFC[common],
                                           "IsProblem" = common %in% problem_genes,
                                           "Tissue" = gsub("(","\n(",tissue,fixed=T)))
}

# Distribution of correlations for young and old
pdf("Plots/5_Predictability/young_old_cordist.pdf", width = 15, height = 10)
ggplot(plot.data) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_density(aes(x = YoungCor, color = "Young")) +
  geom_density(aes(x = OldCor, color = "Old")) +
  scale_color_manual(values = c("Old" = "#FB8500", "Young" = "#126782"))+
  facet_wrap(~ Tissue) + xlab("Correlations (predictability)") +
  labs(color = "Age") + coord_cartesian(xlim = c(-1,1)) +
  theme(text = element_text(size = 20), legend.position = "bottom",
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
dev.off()


# denominator is always negative
denominator <- (log(1-plot.data$OldCor)+log(1-plot.data$YoungCor))/2

png("Plots/5_Predictability/denominator_corLFC.png")
plot(density(denominator), main = "Denominators across all genes / tissues")
dev.off()

adj_denominator <- (log2(2-plot.data$OldCor)+log2(2-plot.data$YoungCor))/2

png("Plots/5_Predictability/denominator_corLFC.png")
plot(density(adj_denominator), main = "Denominators across all genes / tissues")
dev.off()


adj.plot.data <- plot.data
adj.plot.data$PrevLFCs[adj.plot.data$PrevLFCs > 5] <- 5
adj.plot.data$PrevLFCs[adj.plot.data$PrevLFCs < -5] <- -5

pdf("Plots/5_Predictability/young_old_corscatter.pdf", width = 15, height = 10)
ggplot(adj.plot.data) +
  geom_point(aes(x = YoungCor, y = OldCor, color = PrevLFCs),
             size = .1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_color_gradient2(low = "#126782", mid = "white",
                        high = "#FB8500", midpoint = 0) +
  facet_wrap(~ Tissue, scales = "free") +
  labs(color = "Age LFC") +
  xlim(c(-1,1)) + ylim(c(-1,1)) +
  xlab("Correlation in young group") + ylab("Correlation in old group") +
  theme(text = element_text(size = 20), legend.position = "bottom")
dev.off()


pdf("Plots/5_Predictability/young_old_corscatter_fixed.pdf", width = 15, height = 10)
ggplot(adj.plot.data) +
  geom_point(aes(x = YoungCor, y = OldCor, color = AdjLFCs),
             size = .1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_color_gradient2(low = "#126782", mid = "white",
                        high = "#FB8500", midpoint = 0) +
  facet_wrap(~ Tissue) +
  labs(color = "Age LFC") +
  xlim(c(-1,1)) + ylim(c(-1,1)) +
  xlab("Correlation in young group") + ylab("Correlation in old group") +
  theme(text = element_text(size = 20), legend.position = "bottom")
dev.off()

# this one
ggplot(adj.plot.data) +
  geom_point(aes(x = OldCor-YoungCor, y = AdjLFCs,
                 color = (OldCor+YoungCor)/2),
             size = .1) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed") +
  scale_color_viridis_c(direction = -1) +
  xlim(c(-1,1)) + ylim(c(-1,1)) + labs(color = "Average cor.") +
  facet_wrap(~ Tissue)
  
# this one
ggplot(adj.plot.data) +
  geom_point(aes(x = PrevLFCs, y = AdjLFCs,
                 color = OldCor-YoungCor),
             size = .1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_gradient2(high = "#126782", mid = "white",
                        low = "#FB8500", midpoint = 0) +
  # xlim(c(-1,1)) + ylim(c(-1,1)) + labs(color = "Average cor.") +
  facet_wrap(~ Tissue) + labs(color = "Old - Young cor.") +
  theme(text = element_text(size = 20), legend.position = "bottom")



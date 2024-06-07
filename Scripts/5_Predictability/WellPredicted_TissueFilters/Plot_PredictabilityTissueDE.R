
source("Scripts/5_Predictability/params.R")
source("Scripts/5_Predictability/PoorlyPredicted/params.R")
source("Scripts/functions.R")

predictability_files <- list.files("Outputs/5_Predictability/WellPredicted_TissueFilters",
                                   pattern = "sampled_centered_correlations_spearman",
                                   full.names = TRUE)
DE_files <- list.files("Outputs/5_Predictability/WellPredicted_TissueFilters/TissueDE",
                       pattern = "sampled_tissueDE", full.names = TRUE)
cor_thre <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/correlation_thresholds.rds")

plot.data <- data.frame("Predictability" = NULL, "logFC" = NULL, "Gene" = NULL,
                        "WellPredicted" = NULL, "PoorlyPredicted" = NULL,
                        "Tissue" = NULL)

for(t in ALL_TISSUES){
  
  predictability <- readRDS(predictability_files[grep(t, predictability_files,
                                                      fixed = TRUE)])
  fit <- readRDS(DE_files[grep(t, DE_files, fixed = TRUE)])
  DE <- topTable(fit, coef = paste0("Group",t), number = nrow(fit$coefficients))
  DE <- DataENSGToSymbol(as.matrix(DE))
  genes <- intersect(names(LFC), names(predictability))
  
  plot.data <- rbind.data.frame(plot.data,
                                data.frame("Predictability" = predictability[genes],
                                           "logFC" = LFC[genes],
                                           "Gene" = genes,
                                           "WellPredicted" = predictability[genes] > cor_thre[t],
                                           "PoorlyPredicted" = predictability[genes] < COR_THRE,
                                           "Tissue" = t))
}

library(ggplot2)

ggplot(plot.data) +
  geom_point(aes(x = logFC, y = Predictability), size = .1, alpha = .2) +
  facet_wrap(~ Tissue)

pdf("Plots/5_Predictability/WellPredicted_TissueFilters/predictability_tissueDE.pdf",
    width = 12)
ggplot() +
  geom_density(data = subset(plot.data, !PoorlyPredicted),
               aes(x = logFC, color = "Background"), fill = "grey") +
  geom_density(data = subset(plot.data, PoorlyPredicted),
               aes(x = logFC, color = "Poorly predicted genes"), size = 1) +
  xlab("logFC (tissue vs all others)") +
  scale_color_manual(values = c("Poorly predicted genes" = POORLY_PRED_COL,
                                "Background" = "black"),
                     name = NULL) +
  facet_wrap(~ Tissue, labeller = labeller(Tissue = tissue_labels_dash)) +
  theme_classic() + theme(text = element_text(size = 20),
                          axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          legend.position = "bottom")
dev.off()


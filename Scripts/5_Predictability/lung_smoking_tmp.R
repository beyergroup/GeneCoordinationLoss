plot.data <- data.frame("NonSmoker" = readRDS("Outputs/3_GTExDataPrep/Smoking/Lung_nonsmoker_sampled_centered_correlations.rds"),
                        "Smoker" = readRDS("Outputs/3_GTExDataPrep/Smoking/Lung_smoker_sampled_centered_correlations.rds"))

well_predicted <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")

plot.data <- plot.data[intersect(rownames(plot.data),well_predicted$Lung),]

ggplot(plot.data) +
  geom_density(aes(x = NonSmoker, color = "NonSmoker")) +
  geom_density(aes(x = Smoker, color = "Smoker")) +
  scale_color_manual(values = c("NonSmoker" = "darkgreen",
                                "Smoker" = "orangered"),
                     name = "Data subset") + xlab("Pearson correlation") +
  theme_classic() + theme(text = element_text(size = 20),
                          axis.title.y = element_blank())
ggsave("Plots/5_Predictability/tmp/smoking_predictability.png")


ggplot(plot.data) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 1) +
  geom_abline(slope = 1, intercept = 0.2, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -0.2, linetype = "dashed") +
  geom_point(aes(x = NonSmoker, y = Smoker), alpha = .2, size = 1) +
  xlab("Pearson cor. (non-smokers)") + ylab("Pearson cor. (smokers)") +
  theme_classic() + theme(text = element_text(size = 20))
ggsave("Plots/5_Predictability/tmp/smoking_predictability_scatter.png")

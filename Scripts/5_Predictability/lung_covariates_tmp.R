
library(reshape2)
library(ggplot2)

data <- list("NoVentilator" = readRDS("Outputs/3_GTExDataPrep/Ventilator/Lung_noventilator_ageFC_0.75.rds")[,"Slope",drop=F],
             "Ventilator" = readRDS("Outputs/3_GTExDataPrep/Ventilator/Lung_ventilator_ageFC_0.75.rds")[,"Slope",drop=F])
data <- lapply(data, scale, center = F)
plot.data <- melt(data)

global <- as.data.frame(readRDS("Outputs/5_Predictability/tmp/Lung_ageslope_0.75.rds"))
global$Slope <- scale(global$Slope, center = F)

ggplot(plot.data) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_density(data = global,
               aes(x = Slope, color = "Regression data")) +
  geom_density(aes(x = value, color = L1)) + xlab("Slope / Fold Change") +
  scale_color_manual(name = "Data subset",
                     values = c("NoVentilator" = "darkgreen",
                                "Regression data" = "black",
                                "Ventilator" = "orangered")) +
  theme_classic() + theme(text = element_text(size = 20),
                          axis.title.y = element_blank())
ggsave("Plots/5_Predictability/tmp/Lung_ventilator_scaled_ageslopes.png")

plot.data$Global <- global[as.character(plot.data$Var1),"Slope"]

ggplot(plot.data) +
  geom_point(aes(x = Global, y = value)) +
  geom_smooth(aes(x = Global, y = value), method = "lm") +
  xlab("Regression slope") + ylab("Subset FC") +
  facet_wrap(~ L1) + theme(text = element_text(size = 20))
ggsave("Plots/5_Predictability/tmp/Lung_ventilator_ageslopeFC.png")

rm(data,plot.data,global); gc()


data <- list("NonSmoker" = readRDS("Outputs/3_GTExDataPrep/Smoking/Lung_nonsmoker_ageFC_0.75.rds")[,"Slope",drop=F],
             "Smoker" = readRDS("Outputs/3_GTExDataPrep/Smoking/Lung_smoker_ageFC_0.75.rds")[,"Slope",drop=F])
data <- lapply(data, scale, center = F)
plot.data <- melt(data)

global <- as.data.frame(readRDS("Outputs/5_Predictability/tmp/Lung_ageslope_0.75.rds"))
global$Slope <- scale(global$Slope, center = F)

ggplot(plot.data) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_density(data = global,
               aes(x = Slope, color = "Regression data")) +
  geom_density(aes(x = value, color = L1)) +  xlab("Slope / Fold Change") +
  scale_color_manual(name = "Data subset",
                     values = c("NonSmoker" = "darkgreen",
                                "Regression data" = "black",
                                "Smoker" = "orangered")) +
  theme_classic() + theme(text = element_text(size = 20),
                          axis.title.y = element_blank())
ggsave("Plots/5_Predictability/tmp/Lung_smoker_scaled_ageslopes.png")

plot.data$Global <- global[as.character(plot.data$Var1),"Slope"]

ggplot(plot.data) +
  geom_point(aes(x = Global, y = value)) +
  geom_smooth(aes(x = Global, y = value), method = "lm") +
  xlab("Regression slope") + ylab("Subset FC") +
  facet_wrap(~ L1) + theme(text = element_text(size = 20))
ggsave("Plots/5_Predictability/tmp/Lung_smoker_ageslopeFC.png")

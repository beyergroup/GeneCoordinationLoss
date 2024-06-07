
library(ggplot2)
library(ggpubr)
source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")

EXPR_THRE = 0.001
VAR_THRE = 0.005


joint_data <- ReadRDS("Outputs/5_Predictability/Age/Age_MaxSubset/joint_data.rds")

plot.data <- data.frame()

for(tissue in unique(joint_data$Tissue)){
  
  data <- subset(joint_data, Tissue == tissue)
  
  plot.data <- rbind.data.frame(plot.data,
                                data.frame(
                                  "Scenarios" = c(rep("Differential\nexpression",4),
                                                  rep("Differential\nvariance",4),
                                                  rep("No\nchange",2),
                                                  rep("Total",2)),
                                  "Predictability" = rep(c("Up","Down"),6),
                                  "Cases" = c(rep("Decrease",2),
                                              rep("Increase",2),
                                              rep("Decrease",2),
                                              rep("Increase",2),
                                              rep("",4)),
                                  "Count" = c(sum((data$DE_LFC < -EXPR_THRE) & (data$IsTopHit)),
                                              sum((data$DE_LFC < -EXPR_THRE) & (data$IsBottomHit)),
                                              sum((data$DE_LFC > EXPR_THRE) & (data$IsTopHit)),
                                              sum((data$DE_LFC > EXPR_THRE) & (data$IsBottomHit)),
                                      
                                              sum((data$DV_Slope < -VAR_THRE) & (data$IsTopHit)),
                                              sum((data$DV_Slope < -VAR_THRE) & (data$IsBottomHit)),
                                              sum((data$DV_Slope > VAR_THRE) & (data$IsTopHit)),
                                              sum((data$DV_Slope > VAR_THRE) & (data$IsBottomHit)),
                                      
                                              sum((abs(data$DE_LFC) < EXPR_THRE) & (abs(data$DV_Slope) < VAR_THRE) & (data$IsTopHit)),
                                              sum((abs(data$DE_LFC) < EXPR_THRE) & (abs(data$DV_Slope) < VAR_THRE) & (data$IsBottomHit)),
                                              
                                              sum(data$IsTopHit),
                                              sum(data$IsBottomHit)),
                                  "Tissue" = tissue))
  
  
}

p <- ggplot(plot.data) +
  geom_bar(aes(x = Cases, fill = Predictability, y = Count),
           stat = "identity", position = "dodge", width = .6, color = "white") +
  coord_flip() + xlab("") + ylab("Gene number") +
  scale_fill_manual(values =
                      c("Down" = unname(age_palette["60-69"]),
                        "Up" = unname(age_palette["20-29"])),
                    labels = c("Up" = "Increasing predictability w/ age",
                               "Down" = "Decreasing predictability w/ age"),
                    name = NULL) +
  facet_grid(Scenarios ~ Tissue, scales = "free", switch = "y", space = "free",
             labeller = labeller(Tissue = tissue_labels_dash)) +
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.position = "top", strip.background.y = element_blank(),
        strip.placement = "outside", axis.ticks.y = element_blank())


TISSUE = "WholeBlood"


simple_scat <- ggplot(subset(joint_data, Tissue == TISSUE)) +
  # axes
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  # background points
  geom_point(data = subset(joint_data, (Tissue == TISSUE) & (!IsTopHit) &
                             (!IsBottomHit)),
             aes(x = DE_LFC, y = DV_Slope), size = .5, color = "grey") +
  # predictability hits
  geom_point(data = subset(joint_data, (Tissue == TISSUE) & IsTopHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Increasing predictability w/ age")) +
  geom_point(data = subset(joint_data, (Tissue == TISSUE) & IsBottomHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Decreasing predictability w/ age")) +
  # adjustments
  xlab("Expression slope") + ylab("Variance slope") + guides(color = "none") +
  scale_color_manual(values =
                       c("Decreasing predictability w/ age" =
                           unname(age_palette["60-69"]),
                         "Increasing predictability w/ age" =
                           unname(age_palette["20-29"])),
                     name = NULL) +
  theme_classic() + theme(text = element_text(size= 20),
                          plot.title = element_text(hjust = 0.5),
                          axis.text = element_blank(),
                          axis.ticks = element_blank())

simple_dens_expr <- ggplot(subset(joint_data, Tissue == TISSUE)) +
  geom_density(aes(x = DE_LFC), fill = "grey", color = "black") +
  theme_classic() + theme(text = element_text(size = 20),
                          axis.text = element_blank(),
                          axis.ticks = element_blank(),
                          axis.title = element_blank(),
                          axis.line = element_blank())

simple_dens_var <- ggplot(subset(joint_data, Tissue == TISSUE)) +
  geom_density(aes(x = DV_Slope), fill = "grey", color = "black") + coord_flip() +
  theme_classic() + theme(text = element_text(size = 20),
                          axis.text = element_blank(),
                          axis.ticks = element_blank(),
                          axis.title = element_blank(),
                          axis.line = element_blank())

global_scat <- ggarrange(plotlist = list(simple_dens_expr,NULL,
                                         simple_scat,simple_dens_var),
                         heights = c(1,4), widths = c(4,1), align = "hv")

scat_expr <- ggplot(subset(joint_data, Tissue == TISSUE)) +
  # axes
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  # background points
  geom_point(data = subset(joint_data, (Tissue == TISSUE) & (!IsTopHit) &
                             (!IsBottomHit)),
             aes(x = DE_LFC, y = DV_Slope), size = .5, color = "grey") +
  # predictability hits
  geom_point(data = subset(joint_data, (Tissue == TISSUE) & IsTopHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Increasing predictability w/ age")) +
  geom_point(data = subset(joint_data, (Tissue == TISSUE) & IsBottomHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Decreasing predictability w/ age")) +
  # windows for expression thresholds
  geom_rect(aes(xmax = max(DE_LFC), ymin = min(DV_Slope), ymax = max(DV_Slope)),
            xmin = EXPR_THRE, fill = NA, linetype = "dashed", color = "black") +
  geom_rect(aes(xmin = min(DE_LFC), ymin = min(DV_Slope), ymax = max(DV_Slope)),
            xmax = -EXPR_THRE, fill = NA, linetype = "dashed", color = "black") +
  ggtitle("Differential expression") +
  # adjustments
  xlab("Expression slope") + ylab("Variance slope") + guides(color = "none") +
  scale_color_manual(values =
                       c("Decreasing predictability w/ age" =
                           unname(age_palette["60-69"]),
                         "Increasing predictability w/ age" =
                           unname(age_palette["20-29"])),
                     name = NULL) +
  theme_classic() + theme(text = element_text(size= 20),
                          plot.title = element_text(hjust = 0.5),
                          axis.text = element_blank(),
                          axis.ticks = element_blank())


scat_var <- ggplot(subset(joint_data, Tissue == TISSUE)) +
  # axes
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  # background points
  geom_point(data = subset(joint_data, (Tissue == TISSUE) & (!IsTopHit) &
                             (!IsBottomHit)),
             aes(x = DE_LFC, y = DV_Slope), size = .5, color = "grey") +
  # predictability hits
  geom_point(data = subset(joint_data, (Tissue == TISSUE) & IsTopHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Increasing predictability w/ age")) +
  geom_point(data = subset(joint_data, (Tissue == TISSUE) & IsBottomHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Decreasing predictability w/ age")) +
  # windows for variance thresholds
  geom_rect(aes(ymax = max(DV_Slope), xmin = min(DE_LFC), xmax = max(DE_LFC)),
            ymin = VAR_THRE, fill = NA, linetype = "dashed", color = "black") +
  geom_rect(aes(ymin = min(DV_Slope), xmin = min(DE_LFC), xmax = max(DE_LFC)),
            ymax = -VAR_THRE, fill = NA, linetype = "dashed", color = "black") +
  ggtitle("Differential variance") +
  # adjustments
  xlab("Expression slope") + ylab("Variance slope") + guides(color = "none") +
  scale_color_manual(values =
                       c("Decreasing predictability w/ age" =
                           unname(age_palette["60-69"]),
                         "Increasing predictability w/ age" =
                           unname(age_palette["20-29"])),
                     name = NULL) +
  theme_classic() + theme(text = element_text(size= 20),
                          plot.title = element_text(hjust = 0.5),
                          axis.text = element_blank(),
                          axis.ticks = element_blank())

no_scat <- ggplot(subset(joint_data, Tissue == TISSUE)) +
  # axes
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  # background points
  geom_point(data = subset(joint_data, (Tissue == TISSUE) & (!IsTopHit) &
                             (!IsBottomHit)),
             aes(x = DE_LFC, y = DV_Slope), size = .5, color = "grey") +
  # predictability hits
  geom_point(data = subset(joint_data, (Tissue == TISSUE) & IsTopHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Increasing predictability w/ age")) +
  geom_point(data = subset(joint_data, (Tissue == TISSUE) & IsBottomHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Decreasing predictability w/ age")) +
  # windows for both thresholds not met
  geom_rect(ymax = VAR_THRE, xmin = -EXPR_THRE, xmax = EXPR_THRE, ymin = -VAR_THRE,
            fill = NA, linetype = "dashed", color = "black") +
  ggtitle("No change") +
  # adjustments
  xlab("Expression slope") + ylab("Variance slope") + guides(color = "none") +
  scale_color_manual(values =
                       c("Decreasing predictability w/ age" =
                           unname(age_palette["60-69"]),
                         "Increasing predictability w/ age" =
                           unname(age_palette["20-29"])),
                     name = NULL) +
  theme_classic() + theme(text = element_text(size= 20),
                          plot.title = element_text(hjust = 0.5),
                          axis.ticks = element_blank(),
                          axis.text = element_blank())

pdf("Plots/Figures/Parts/Figure3.pdf", width = 21, height = 10)
A <- ggarrange(plotlist = list(global_scat,NULL,scat_expr,scat_var,no_scat),
               nrow = 1, align = "h")
ggarrange(plotlist = list(A,p), labels = c("A","B"), nrow = 2, ncol = 1,
          heights = c(2,3), font.label = list(size = 22))
dev.off()

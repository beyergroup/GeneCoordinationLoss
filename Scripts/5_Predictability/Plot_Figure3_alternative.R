
library(ggplot2)
library(ggpubr)
source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")

SLOPE_DIR = "Outputs/5_Predictability/Age/Age_MaxSubset"
TISSUES = "main"
EXPR_THRE = 0.001
VAR_THRE = 0.001
INCL_70 = F

# Scatter plot expression variance (panel E)

if(INCL_70){
  plot.data <- ReadRDS(paste0(SLOPE_DIR,"/joint_data.rds"))
} else{
  plot.data <- ReadRDS(paste0(SLOPE_DIR,"/joint_data_no70.rds"))
}
if(TISSUES == "all"){
  plot.data$Tissue <- factor(plot.data$Tissue,
                             levels = unique(plot.data$Tissue)[c(1:3,16,4:15)])
} else{
  plot.data$Tissue <- factor(as.character(plot.data$Tissue),
                             levels = MAIN_TISSUES)
}

# plot.data$DV_Slope[which.min(plot.data$DV_Slope)] <- NA
# plot.data$DV_Slope[which.min(plot.data$DV_Slope)] <- NA

RemoveOutliers <- function(df,x){
  
  aux <- df[,x]
  posmed <- median(aux[aux >= 0])
  posmad <- mad(aux[aux >= 0])
  negmed <- median(aux[aux < 0])
  negmad <- mad(aux[aux < 0])
  
  ind <- (df[,x] > (posmed + (20*posmad))) |
    (df[,x] < (negmed - (20*negmad)))
  
  return(df[!ind,])
}

tmp <- do.call(rbind.data.frame,
               sapply(unique(plot.data$Tissue),
                      function(t) RemoveOutliers(subset(plot.data, Tissue == t),
                                                 "DV_Slope"),
                      simplify = F))

scat <- ggplot(tmp) +
  geom_point(data = subset(tmp, (!IsTopHit) & (!IsBottomHit)),
             aes(x = DE_LFC, y = DV_Slope), size = .5, color = "grey") +
  geom_point(data = subset(tmp, IsTopHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Increasing predictability w/ age")) +
  geom_point(data = subset(tmp, IsBottomHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Decreasing predictability w/ age")) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  xlab("Age-related expression changes") +
  ylab("Age-related variance changes") +
  facet_wrap(~ Tissue, nrow = 2, scales = "free",
             labeller = labeller(Tissue = tissue_labels_dash)) +
  scale_color_manual(values =
                       c("Decreasing predictability w/ age" =
                           unname(age_palette["60-69"]),
                         "Increasing predictability w/ age" =
                           unname(age_palette["20-29"])),
                     name = NULL) +
  theme_classic() + theme(text = element_text(size = 20),
                          legend.position = "bottom",
                          axis.text.x = element_text(angle = 40, vjust = 1,
                                                     hjust = 1))

scat.full <- ggplot(plot.data) +
  geom_point(data = subset(plot.data, (!IsTopHit) & (!IsBottomHit)),
             aes(x = DE_LFC, y = DV_Slope), size = .5, color = "grey") +
  geom_point(data = subset(plot.data, IsTopHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Increasing predictability w/ age")) +
  geom_point(data = subset(plot.data, IsBottomHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Decreasing predictability w/ age")) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  xlab("Age-related expression changes") +
  ylab("Age-related variance changes") +
  facet_wrap(~ Tissue, nrow = 2, scales = "free",
             labeller = labeller(Tissue = tissue_labels_dash)) +
  scale_color_manual(values =
                       c("Decreasing predictability w/ age" =
                           unname(age_palette["60-69"]),
                         "Increasing predictability w/ age" =
                           unname(age_palette["20-29"])),
                     name = NULL) +
  theme_classic() + theme(text = element_text(size = 20),
                          legend.position = "bottom",
                          axis.text.x = element_text(angle = 40, vjust = 1,
                                                     hjust = 1))

if(INCL_70){
  WriteRDS(scat, "Outputs/5_Predictability/Age/Age_MaxSubset/scatterplot_crop_main.rds")
  WriteRDS(scat.full, "Outputs/5_Predictability/Age/Age_MaxSubset/scatterplot_main.rds")
  pdf("Plots/Figures/variance_expression_full.pdf", width = 14, height = 6)
  scat.full
  dev.off()
} else{
  WriteRDS(scat, "Outputs/5_Predictability/Age/Age_MaxSubset/scatterplot_crop_main_no70.rds")
  WriteRDS(scat.full, "Outputs/5_Predictability/Age/Age_MaxSubset/scatterplot_main_no70.rds")
  pdf("Plots/Figures/variance_expression_no70.pdf", width = 14, height = 6)
  scat.full
  dev.off()
}

rm(plot.data); gc()



joint_data <- ReadRDS("Outputs/5_Predictability/Age/Age_MaxSubset/joint_data.rds")
joint_data_no70 <- ReadRDS("Outputs/5_Predictability/Age/Age_MaxSubset/joint_data_no70.rds")

plot.data <- data.frame()

for(tissue in unique(joint_data$Tissue)){
  
  data <- subset(joint_data, Tissue == tissue)
  data_no70 <- subset(joint_data_no70, Tissue == tissue)
  
  # subset to common
  genes <- intersect(data$Gene,data_no70$Gene)
  data <- data[match(genes,data$Gene),]
  data_no70 <- data_no70[match(genes,data_no70$Gene),]
  
  hits <- ReadRDS(paste0("Outputs/5_Predictability/Age/Age_MaxSubset/",
                         tissue,"_predictability_hits_robust.rds"))
  
  pred_cor <- ReadRDS(paste0("Outputs/5_Predictability/Age/Age_MaxSubset/",
                             tissue,"_predcor_robust_slope.rds"))
  
  plot.data <- rbind.data.frame(plot.data,
                                data.frame(
                                  "Scenarios" = c(rep("Expression\nchange",4),
                                                  rep("Variance\nchange",4),
                                                  rep("Correlation\nchange",4),
                                                  rep("Total",2)),
                                  "Predictability" = rep(c("Up","Down"),7),
                                  "Cases" = c(rep(c("Decrease","Decrease",
                                                    "Increase","Increase"),3),
                                              rep("",2)),
                                  "Count" = c(sum((data$DE_LFC < -EXPR_THRE) & (data_no70$DE_LFC < -EXPR_THRE) & (data$Gene %in% hits$Up)),
                                              sum((data$DE_LFC < -EXPR_THRE) & (data_no70$DE_LFC < -EXPR_THRE) & (data$Gene %in% hits$Down)),
                                              sum((data$DE_LFC > EXPR_THRE) & (data_no70$DE_LFC > EXPR_THRE) & (data$Gene %in% hits$Up)),
                                              sum((data$DE_LFC > EXPR_THRE) & (data_no70$DE_LFC > EXPR_THRE) & (data$Gene %in% hits$Down)),
                                      
                                              sum((data$DV_Slope < -VAR_THRE) & (data_no70$DV_Slope < VAR_THRE) & (data$Gene %in% hits$Up)),
                                              sum((data$DV_Slope < -VAR_THRE) & (data_no70$DV_Slope < VAR_THRE) & (data$Gene %in% hits$Down)),
                                              sum((data$DV_Slope > VAR_THRE) & (data_no70$DV_Slope > VAR_THRE)  & (data$Gene %in% hits$Up)),
                                              sum((data$DV_Slope > VAR_THRE) & (data_no70$DV_Slope > VAR_THRE) & (data$Gene %in% hits$Down)),
                                      
                                              sum((pred_cor$AvgSlope[match(data$Gene, pred_cor$Gene)] < 0) & (data$Gene %in% hits$Up), na.rm = T),
                                              sum((pred_cor$AvgSlope[match(data$Gene, pred_cor$Gene)] < 0) & (data$Gene %in% hits$Down), na.rm = T),
                                              sum((pred_cor$AvgSlope[match(data$Gene, pred_cor$Gene)] > 0) & (data$Gene %in% hits$Up), na.rm = T),
                                              sum((pred_cor$AvgSlope[match(data$Gene, pred_cor$Gene)] > 0) & (data$Gene %in% hits$Down), na.rm = T),
                                              
                                              sum(data$IsTopHit),
                                              sum(data$IsBottomHit)),
                                  "Tissue" = tissue))
}

plot.data$Scenarios <- factor(as.character(plot.data$Scenarios),
                              levels = c("Expression\nchange",
                                         "Variance\nchange",
                                         "Correlation\nchange",
                                         "Total"))

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


pdf("Plots/Figures/Parts/Figure3_B_crop.pdf", width = 14, height = 6)
# A <- ggarrange(plotlist = list(global_scat,NULL,scat_expr,scat_var,NULL,NULL),
#                nrow = 1, align = "h")
ggarrange(scat, labels = "B", font.label = list(size = 22))
dev.off()

pdf("Plots/Figures/Parts/Figure3_C.pdf", width = 20, height = 6)
ggarrange(p, labels = "C", font.label = list(size = 22))
dev.off()

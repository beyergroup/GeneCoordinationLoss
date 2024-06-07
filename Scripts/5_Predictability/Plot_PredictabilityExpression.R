# Plot delta LFCs vs expression LFCs

args = commandArgs(trailingOnly=TRUE)
COR_DIR = "Outputs/5_Predictability/Age/Age_MaxSubset"
DE_DIR = "Outputs/3_GTExDataPrep/Differential_Expression/Max_Subset"
PLOT_DIR = "Plots/5_Predictability/Age/Age_MaxSubset"
PATTERN = "well_predicted"
INCL_70 = F

.libPaths("Resources/Rlibs/R-4.0.3/")
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(viridis)
source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")

if(INCL_70){
  cor_files <- list.files(COR_DIR, pattern = paste0("ageslope_",PATTERN,".rds"),
                          full.names = T)
  expression_files <- list.files(DE_DIR, pattern = "_DE.rds", full.names = T)
  var_files <- list.files(DE_DIR, pattern = "_var_regression.rds", full.names = T)
} else{
  cor_files <- list.files(COR_DIR, pattern = paste0("ageslope_",
                                                    PATTERN,"_no70.rds"),
                          full.names = T)
  expression_files <- list.files(DE_DIR, pattern = "_DE_no70.rds", full.names = T)
  var_files <- list.files(DE_DIR, pattern = "_var_regression_no70.rds",
                          full.names = T)
}

plots <- list()
plot.data <- data.frame()

for(tissue in MAIN_TISSUES){
  
  ageDE <- ReadRDS(expression_files[grepl(tissue,expression_files,fixed = T)])
  ageDE <- limma::topTable(fit = ageDE, coef = "Age", number = nrow(ageDE))
  ageDE <- DataENSGToSymbol(ageDE, remove_dup = T)
  
  ageDV <- ReadRDS(var_files[grepl(tissue,var_files,fixed=T)])
  ageDV <- DataENSGToSymbol(ageDV, remove_dup = T)
  
  ageDP <- ReadRDS(cor_files[grepl(tissue,cor_files,fixed = T)])
  
  common <- intersect(intersect(rownames(ageDP),rownames(ageDE)),
                      rownames(ageDV))
  
  message(tissue)
  message("mean = ",round(mean(ageDP[common,"Slope"]),2))
  message("median = ",round(median(ageDP[common,"Slope"]),2),"\n")
  
  if(INCL_70){
    hits <- readRDS(paste0(COR_DIR,"/",tissue,"_predictability_hits.rds"))
  } else{
    hits <- readRDS(paste0(COR_DIR,"/",tissue,"_predictability_hits_no70.rds"))
  }
  
  curr.plot.data <- data.frame("Gene" = common,
                               "DE_LFC" = ageDE[common,"logFC"],
                               "DP_Slope" = ageDP[common,"Slope"],
                               "DE_pval" = ageDE[common,"P.Value"],
                               "DV_Slope" = ageDV[common,"Slope"],
                               "DV_pval" = ageDV[common,"pval"],
                               "Tissue" = tissue,
                               "IsTopHit" = common %in% hits$Up,
                               "IsBottomHit" = common %in% hits$Down)
  plot.data <- rbind.data.frame(plot.data, curr.plot.data)
  
  plots[[tissue]] <- ggplot(curr.plot.data) +
    geom_density(aes(x = DV_Slope, fill = "Background"), color = "transparent") +
    geom_density(data = subset(curr.plot.data, IsTopHit),
                 aes(x = DV_Slope, color = "Increasing correlation w/ age"),
                 size = 1) +
    geom_density(data = subset(curr.plot.data, IsBottomHit),
                 aes(x = DV_Slope, color = "Decreasing correlation w/ age"),
                 size = 1) +
    xlab("Variance slope") + ggtitle(tissue_labels_dash[tissue]) +
    scale_color_manual(values =
                         c("Decreasing correlation w/ age" =
                             unname(age_palette["60-69"]),
                           "Increasing correlation w/ age" =
                             unname(age_palette["20-29"])),
                       name = NULL) +
    scale_fill_manual(values = c("Background" = "grey"), name = NULL) +
    theme_classic() + theme(text = element_text(size = 20),
                            legend.position = "bottom",
                            axis.title.y = element_blank(),
                            axis.text.y = element_blank(),
                            axis.ticks.y = element_blank())
}

if(INCL_70){
  WriteRDS(plot.data,
           "Outputs/5_Predictability/Age/Age_MaxSubset/joint_data.rds")
} else{
  WriteRDS(plot.data,
           "Outputs/5_Predictability/Age/Age_MaxSubset/joint_data_no70.rds")
}


if(INCL_70){
  pdf(paste0(PLOT_DIR,"/slopes_cor_expr_",PATTERN,".pdf"),
      width = 12, height = 6)
} else{
  pdf(paste0(PLOT_DIR,"/slopes_cor_expr_",PATTERN,"_no70.pdf"),
      width = 12, height = 6)
}
ggplot(plot.data) +
  geom_point(data = subset(plot.data, (DE_pval >= 0.05) &
                             (!IsTopHit) & (!IsBottomHit)),
             aes(x = DE_LFC, y = DP_Slope), color = "grey", size = .5) +
  geom_point(data = subset(plot.data, DE_pval < 0.05),
             aes(x = DE_LFC, y = DP_Slope), color = "darkgrey") +
  geom_point(data = subset(plot.data, IsBottomHit),
             aes(x = DE_LFC, y = DP_Slope,
                 color = "Decreasing correlation w/ age")) +
  geom_point(data = subset(plot.data, IsTopHit),
             aes(x = DE_LFC, y = DP_Slope,
                 color = "Increasing correlation w/ age")) +
  scale_color_manual(values =
                       c("Decreasing correlation w/ age" =
                           unname(age_palette["60-69"]),
                         "Increasing correlation w/ age" =
                           unname(age_palette["20-29"])),
                     name = NULL) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  xlim(c(-max(abs(plot.data$DE_LFC)), max(abs(plot.data$DE_LFC)))) +
  xlab("Expression slope") + ylab("Correlation slope") +
  facet_wrap(~ Tissue, nrow = 2, labeller = labeller(Tissue = tissue_labels_dash)) +
  theme_classic() +
  theme(text = element_text(size = 20), title = element_text(size = 18),
        legend.position = "bottom")
dev.off()

if(INCL_70){
  pdf(paste0(PLOT_DIR,"/slopes_cor_var_",PATTERN,".pdf"),
      width = 12, height = 6)
} else{
  pdf(paste0(PLOT_DIR,"/slopes_cor_var_",PATTERN,"_no70.pdf"),
      width = 12, height = 6)
}
ggplot(plot.data) +
  geom_point(data = subset(plot.data, (DV_pval >= 0.005) &
                             (!IsTopHit) & (!IsBottomHit)),
             aes(x = DV_Slope, y = DP_Slope), color = "grey", size = .5) +
  geom_point(data = subset(plot.data, DV_pval < 0.005),
             aes(x = DV_Slope, y = DP_Slope), color = "darkgrey") +
  geom_point(data = subset(plot.data, IsBottomHit),
             aes(x = DV_Slope, y = DP_Slope,
                 color = "Decreasing correlation w/ age")) +
  geom_point(data = subset(plot.data, IsTopHit),
             aes(x = DV_Slope, y = DP_Slope,
                 color = "Increasing correlation w/ age")) +
  scale_color_manual(values =
                       c("Decreasing correlation w/ age" =
                           unname(age_palette["60-69"]),
                         "Increasing correlation w/ age" =
                           unname(age_palette["20-29"])),
                     name = NULL) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  xlim(c(-max(abs(plot.data$DV_Slope)), max(abs(plot.data$DV_Slope)))) +
  xlab("Variance slope") +
  ylab("Correlation slope") +
  facet_wrap(~ Tissue, nrow = 2, labeller = labeller(Tissue = tissue_labels_dash)) +
  theme_classic() +
  theme(text = element_text(size = 20), title = element_text(size = 18),
        legend.position = "bottom")
dev.off()



# Density plots
if(INCL_70){
  pdf(paste0(PLOT_DIR,"/var_densityplots_",PATTERN,".pdf"),
      width = 16, height = 6)
} else{
  pdf(paste0(PLOT_DIR,"/var_densityplots_",PATTERN,"_no70.pdf"),
      width = 16, height = 6)
}
ggarrange(plotlist = plots, common.legend = T, nrow = 2, ncol = 4, legend = "bottom")
dev.off()


if(INCL_70){
  pdf(paste0(PLOT_DIR,"/slopes_scatter_",PATTERN,".pdf"),
      width = 14, height = 6)
} else{
  pdf(paste0(PLOT_DIR,"/slopes_scatter_",PATTERN,"_no70.pdf"),
      width = 14, height = 6)
}
# ggplot(plot.data) +
#   geom_point(aes(x = DE_LFC, y = DV_Slope, color = DP_Slope)) +
#   geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
#   xlab("Expression slope") + ylab("Variance slope") +
#   facet_wrap(~ Tissue, nrow = 2, scales = "free",
#              labeller = labeller(Tissue = tissue_labels_dash)) +
#   scale_color_viridis_c(direction = -1, option = "cividis",
#                         name = "Correlation\nslope") +
#   theme_classic() + theme(text = element_text(size = 20),
#                           axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
ggplot(plot.data) +
  geom_point(data = subset(plot.data, (!IsTopHit) & (!IsBottomHit)),
             aes(x = DE_LFC, y = DV_Slope), size = .5, color = "grey") +
  geom_point(data = subset(plot.data, IsTopHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Increasing correlation w/ age")) +
  geom_point(data = subset(plot.data, IsBottomHit),
             aes(x = DE_LFC, y = DV_Slope,
                 color = "Decreasing correlation w/ age")) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  xlab("Expression slope") + ylab("Variance slope") +
  facet_wrap(~ Tissue, nrow = 2, scales = "free",
             labeller = labeller(Tissue = tissue_labels_dash)) +
  scale_color_manual(values =
                       c("Decreasing correlation w/ age" =
                           unname(age_palette["60-69"]),
                         "Increasing correlation w/ age" =
                           unname(age_palette["20-29"])),
                     name = NULL) +
  theme_classic() + theme(text = element_text(size = 20),
                          legend.position = "bottom",
                          axis.text.x = element_text(angle = 60, vjust = 1,
                                                     hjust = 1))
dev.off()


if(INCL_70){
  pdf(paste0(PLOT_DIR,"/var_volcano_",PATTERN,".pdf"),
      width = 14, height = 6)
} else{
  pdf(paste0(PLOT_DIR,"/var_volcano_",PATTERN,"_no70.pdf"),
      width = 14, height = 6)
}

ggplot(plot.data) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_point(data = subset(plot.data, (!IsTopHit) & (!IsBottomHit)),
             aes(x = DV_Slope, y = -log10(DV_pval)),
             color = "grey", size = .5) +
  geom_point(data = subset(plot.data, IsTopHit),
             aes(x = DV_Slope, y = -log10(DV_pval),
                 color = "Increasing correlation w/ age")) +
  facet_wrap(~ Tissue, nrow = 2, scales = "free",
             labeller = labeller(Tissue = tissue_labels_dash)) +
  xlab("Variance slope") + ylab("Variance -log10(p-value)") +
  scale_color_manual(values =
                       c("Increasing correlation w/ age" =
                           unname(age_palette["20-29"])),
                     name = NULL) +
  theme_classic() + theme(text = element_text(size = 20),
                          legend.position = "bottom")

ggplot(plot.data) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_point(data = subset(plot.data, (!IsTopHit) & (!IsBottomHit)),
             aes(x = DV_Slope, y = -log10(DV_pval)),
             color = "grey", size = .5) +
  geom_point(data = subset(plot.data, IsBottomHit),
             aes(x = DV_Slope, y = -log10(DV_pval),
                 color = "Decreasing correlation w/ age")) +
  facet_wrap(~ Tissue, nrow = 2, scales = "free",
             labeller = labeller(Tissue = tissue_labels_dash)) +
  xlab("Variance slope") + ylab("Variance -log10(p-value)") +
  scale_color_manual(values =
                       c("Decreasing correlation w/ age" =
                           unname(age_palette["60-69"])),
                     name = NULL) +
  theme_classic() + theme(text = element_text(size = 20),
                          legend.position = "bottom")

dev.off()

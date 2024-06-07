# Plot delta LFCs vs expression LFCs

args = commandArgs(trailingOnly=TRUE)
COR_DIR = "Outputs/5_Predictability/Age"
DE_DIR = "Outputs/3_GTExDataPrep/Differential_Expression"
PLOT_DIR = "Plots/5_Predictability/Age"
PATTERN = "well_predicted"
INCL_70 = T

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
  var_files <- list.files(DE_DIR, pattern = "_var_regression", full.names = T)
} else{
  cor_files <- list.files(COR_DIR, pattern = paste0("ageslope_",
                                                    PATTERN,"_no70.rds"),
                          full.names = T)
  expression_files <- list.files(DE_DIR, pattern = "_DE_no70.rds", full.names = T)
  var_files <- list.files(DE_DIR, pattern = "_var_regression_no70.rds",
                          full.names = T)
}

tissues <- unique(sapply(cor_files,
                         function(x) strsplit(tail(strsplit(x,"/")[[1]],1),
                                              "_")[[1]][1]))

plots <- list()
densities <- list()
plot.data <- data.frame()

for(tissue in tissues){
  
  ageDE <- ReadRDS(expression_files[grepl(tissue,expression_files,fixed = T)])
  ageDE <- limma::topTable(fit = ageDE, coef = "AGE", number = nrow(ageDE))
  ageDE <- DataENSGToSymbol(ageDE, remove_dup = T)
  
  ageDV <- ReadRDS(var_files[grepl(tissue,var_files,fixed=T)])
  ageDV <- DataENSGToSymbol(ageDV, remove_dup = T)
  
  ageDP <- ReadRDS(cor_files[grepl(tissue,cor_files,fixed = T)])
  
  common <- intersect(intersect(rownames(ageDP),rownames(ageDE)),
                      rownames(ageDV))
  
  message(tissue)
  message("mean = ",round(mean(ageDP[common,"Slope"]),2))
  message("median = ",round(median(ageDP[common,"Slope"]),2),"\n")
  
  hits <- readRDS(paste0("Outputs/5_Predictability/Age/",
                         tissue,"_predictability_hits.rds"))
  
  curr.plot.data <- data.frame("Gene" = common,
                               "DE_LFC" = ageDE[common,"logFC"],
                               "DP_Slope" = ageDP[common,"Slope"],
                               "DE_pval" = ageDE[common,"adj.P.Val"],
                               "DV_Slope" = ageDV[common,"Slope"],
                               "DV_pval" = ageDV[common,"pval"],
                               "Tissue" = tissue,
                               "IsTopHit" = common %in% hits$Up,
                               "IsBottomHit" = common %in% hits$Down)
  plot.data <- rbind.data.frame(plot.data, curr.plot.data)
  
  d.data <- data.frame("x" = density(curr.plot.data$DP_Slope)$x,
                       "y" = density(curr.plot.data$DP_Slope)$y)
  d.data$Group <- "Center"
  d.data$Group[d.data$x <= -1] <- "LeftTail"
  d.data$Group[d.data$x >= 1] <- "RightTail"
  
  densities[[tissue]] <- ggplot() +
    geom_area(data = d.data,
              aes(x = x, y = y, group = Group, fill = Group),
              color = "black", outline.type = "full") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    annotate("text", x = mean(c(1,max(d.data$x))),
             y = max(d.data$y)/6, col = "#FB8500", size = 5,
             label = paste0("N = ",sum(curr.plot.data$DP_LFC >= 1))) +
    annotate("text", x = mean(c(-1,min(d.data$x))),
             y = max(d.data$y)/6, col = "#126782", size = 5,
             label = paste0("N = ",sum(curr.plot.data$DP_LFC <= -1))) +
    xlab("Error fold change") +
    scale_fill_manual(values = c("LeftTail" = "#126782",
                                 "Center" = "transparent",
                                 "RightTail" = "#FB8500")) +
    guides(fill = "none", color = "none") +
    theme_classic() +
    ggtitle(gsub("(","\n(",tissue,fixed=T)) + 
    theme(text = element_text(size = 20),
          title = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
}

if(INCL_70){
  WriteRDS(plot.data, paste0(COR_DIR,"/joint_data.rds"))
} else{
  WriteRDS(plot.data, paste0(COR_DIR,"/joint_data_no70.rds"))
}


pdf(paste0(PLOT_DIR,"/slopes_cor_expr_",PATTERN,".pdf"),
    width = 12, height = 8)
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
  xlab("Expression slope") +
  ylab("Correlation slope") +
  facet_wrap(~ Tissue) + theme_classic() +
  theme(text = element_text(size = 20), title = element_text(size = 18),
        legend.position = "bottom")
dev.off()

pdf(paste0(PLOT_DIR,"/slopes_cor_var_",PATTERN,".pdf"),
    width = 12, height = 8)
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
  facet_wrap(~ Tissue) + theme_classic() +
  theme(text = element_text(size = 20), title = element_text(size = 18),
        legend.position = "bottom")
dev.off()

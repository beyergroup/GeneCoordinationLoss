PLOT_DIR = "Plots/5_Predictability/Age/Age_MaxSubset"


library(reshape2)
library(ggplot2)
library(viridis)
source("Scripts/functions.R")

expr_files <- list.files("Outputs/3_GTExDataPrep/Subset_Data/Max_Subset",
                         pattern = "sampled_centered_data.rds",
                         full.names = T)
pred_files <- list.files("Outputs/5_Predictability/Age/Age_MaxSubset",
                         pattern = "sampled_centered_net_predictions.rds",
                         full.names = T)
cor_files <- list.files("Outputs/5_Predictability/Age/Age_MaxSubset",
                        pattern = "sampled_centered_correlations_spearman.rds",
                        full.names = T)

metadata <- ReadRDS("Outputs/3_GTExDataPrep/metadata_restrictedaccess.rds")

joint_data <- ReadRDS("Outputs/5_Predictability/Age/Age_MaxSubset/joint_data.rds")
joint_data_no70 <- ReadRDS("Outputs/5_Predictability/Age/Age_MaxSubset/joint_data_no70.rds")
joint_data_small <- ReadRDS("Outputs/5_Predictability/Age/joint_data.rds")

for(tissue in unique(joint_data$Tissue)){
  message(tissue)
  common_top_hits <- intersect(joint_data$Gene[(joint_data$IsTopHit) &
                                                 (joint_data$Tissue == tissue)],
                               intersect(joint_data_no70$Gene[(joint_data_no70$IsTopHit) &
                                                                (joint_data_no70$Tissue == tissue)],
                                         joint_data_small$Gene[(joint_data_small$IsTopHit) &
                                                                 (joint_data_small$Tissue == tissue)]))
  
  message("Main data subset:")
  print(subset(joint_data, (Tissue == tissue) & (Gene %in% common_top_hits)))
  message("No 70s subset:")
  print(subset(joint_data_no70, (Tissue == tissue) & (Gene %in% common_top_hits)))
  message("Small data subset:")
  print(subset(joint_data_small, (Tissue == tissue) & (Gene %in% common_top_hits)))
  
  common_bottom_hits <- intersect(joint_data$Gene[(joint_data$IsBottomHit) &
                                                    (joint_data$Tissue == tissue)],
                                  intersect(joint_data_no70$Gene[(joint_data_no70$IsBottomHit) &
                                                                   (joint_data_no70$Tissue == tissue)],
                                            joint_data_small$Gene[(joint_data_small$IsBottomHit) &
                                                                    (joint_data_small$Tissue == tissue)]))
  
  message("Main data subset:")
  print(subset(joint_data, (Tissue == tissue) & Gene %in% common_bottom_hits))
  message("No 70s subset:")
  print(subset(joint_data_no70, (Tissue == tissue) & Gene %in% common_bottom_hits))
  message("Small data subset:")
  print(subset(joint_data_small, (Tissue == tissue) & Gene %in% common_bottom_hits))
}


data <- joint_data
# data <- subset(joint_data, (IsBottomHit) & (DV_pval < 0.01))
# data <- subset(joint_data, (IsTopHit) & (DV_pval < 0.01))
# data <- subset(joint_data, abs(DP_Slope) >= 0.015)

if(INCL_70){
  # pdf(paste0(PLOT_DIR,"/bottom_hits_var_0.01.pdf"), width = 8)
  pdf(paste0(PLOT_DIR,"/top_hits_var_0.01.pdf"), width = 8)
  # pdf(paste0(PLOT_DIR,"/hits.pdf"), width = 8)
} else{
  # pdf(paste0(PLOT_DIR,"/bottom_hits_var_0.01_no70.pdf"), width = 8)
  pdf(paste0(PLOT_DIR,"/top_hits_var_0.01_no70.pdf"), width = 8)
  # pdf(paste0(PLOT_DIR,"/hits_no70.pdf"), width = 8)
}

for(tissue in unique(data$Tissue)){

  # Batch-corrected, centered expression
  tmp_files <- expr_files[grep(tissue,expr_files,fixed = T)]
  centered_expr <- sapply(tmp_files, ReadRDS)
  centered_expr <- lapply(centered_expr, DataENSGToSymbol, remove_dup = T)
  centered_expr <- do.call(cbind, centered_expr)
  tmp_files <- pred_files[grep(tissue,pred_files,fixed = T)]
  centered_pred <- sapply(tmp_files, ReadRDS)
  centered_pred <- do.call(cbind, centered_pred)
  tmp_files <- cor_files[grep(tissue,cor_files,fixed = T)]
  correlations <- sapply(tmp_files, ReadRDS)
  colnames(correlations) <- sapply(colnames(correlations),
                                   function(x) strsplit(tail(strsplit(x,"/")[[1]],
                                                             1), "_")[[1]][2])
  rm(tmp_files)
  
  # Metadata
  m <- metadata[match(colnames(centered_expr),metadata$SAMPID),]
  
  for(gene in subset(data, Tissue == tissue)$Gene){
    plot.data <- data.frame("Expr" = centered_expr[gene,],
                            "Pred" = centered_pred[gene,],
                            "AgeGroup" = m$AGE_GROUP,
                            "Sex" = m$SEX)
    labs <- paste0(colnames(correlations),"\nrho = ",
                   round(correlations[gene,],2))
    names(labs) <- colnames(correlations)

    scat <- ggplot(plot.data) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
            geom_point(aes(x = Expr, y = Pred, color = AgeGroup)) +
            geom_smooth(aes(x = Expr, y = Pred, color = AgeGroup),
                        method = "lm", alpha = .1) +
            scale_color_viridis(discrete = T, option = "D",
                                end = .9, name = NULL) +
            facet_grid(~ AgeGroup,
                       labeller = labeller(AgeGroup = labs)) +
            xlab(paste(gene,"observed expr.")) +
            ylab(paste(gene,"predicted expr.")) +
            theme_classic() + theme(text = element_text(size = 20),
                                    legend.position = "bottom")
    
    WriteRDS(scat, paste0("Outputs/5_Predictability/Age/Age_MaxSubset/",
                          gene,"_",tissue,"_scatterplot_predobs.rds"))
    
    print(scat +
            ggtitle(tissue, subtitle = paste0("Cor. slope = ",
                                              round(data[gene,"DP_Slope"],3),
                                              "\nVar. slope = ",
                                              round(data[gene,"DV_Slope"],3),
                                              "; p-value = ",
                                              round(data[gene,"DV_pval"],4),
                                              "\nExpr. LFC = ",
                                              round(data[gene,"DE_LFC"],3),
                                              "; p-value = ",
                                              round(data[gene,"DE_pval"],4))))

    plot.data <- melt(plot.data, id.vars = c("AgeGroup","Sex"))
    
    print(ggplot(plot.data) +
            geom_boxplot(aes(x = AgeGroup, y = value, fill = variable),
                         outlier.colour = "transparent", outlier.size = 0) +
            geom_point(aes(x = AgeGroup, y = value, color = variable,
                           shape = Sex),
                       position = position_jitterdodge(jitter.width = 0.05)) +
            xlab("Age group") + ylab(paste(gene,"centered expr.")) +
            scale_fill_manual(values = c("Expr" = obs_color,
                                         "Pred" = net_color),
                              name = NULL,
                              labels = c("Expr" = "Observed expr.",
                                         "Pred" = "Predicted expr.")) +
            scale_color_manual(values = c("Expr" = "black",
                                          "Pred" = "black"),
                               name = NULL) + guides(color = "none") +
            scale_shape_manual(values = c("Female" = 16, "Male" = 3),
                               name = NULL) +
            ggtitle(tissue, subtitle = paste0("Cor. slope = ",
                                              round(data[gene,"DP_Slope"],3),
                                              "\nVar. slope = ",
                                              round(data[gene,"DV_Slope"],3),
                                              "; p-value = ",
                                              round(data[gene,"DV_pval"],4),
                                              "\nExpr. LFC = ",
                                              round(data[gene,"DE_LFC"],3),
                                              "; p-value = ",
                                              round(data[gene,"DE_pval"],4))) +
            theme_classic() + theme(text = element_text(size = 20),
                                    legend.position = "bottom"))
  }
}

dev.off()


data <- subset(joint_data, (Tissue == "Adipose-Subcutaneous") &
                 (IsTopHit | IsBottomHit) & (DE_LFC > 0))
data <- head(data[order(data$DV_Slope, decreasing = T),],6)


pdf(paste0(PLOT_DIR,"/hits_",tissue,".pdf"))

# Look at top hits

source("Scripts/functions.R")
library(reshape2)
library(ggplot2)
library(ggpubr)

# MODE = "YvsO_adj"
MODE = "YvsO_adj_spearman"

expr_files <- list.files("Outputs/5_Predictability",
                         pattern = "sampled_centered_data",
                         full.names = T)
young_expr_files <- expr_files[grep("young",expr_files)]
old_expr_files <- expr_files[grep("old",expr_files)]
rm(expr_files); gc()

pred_files <- list.files("Outputs/5_Predictability",
                         pattern = "sampled_centered_net_predictions.rds",
                         full.names = T)
young_pred_files <- pred_files[grep("young",pred_files)]
old_pred_files <- pred_files[grep("old",pred_files)]
rm(pred_files); gc()

cor_files <- list.files("Outputs/5_Predictability",
                        pattern = c("YvsO_adj" = "sampled_centered_correlations.rds",
                                    "YvsO_adj_spearman" = "sampled_centered_correlations_spearman.rds")[MODE],
                        full.names = T)
young_cor_files <- cor_files[grep("young",cor_files)]
old_cor_files <- cor_files[grep("old",cor_files)]
rm(cor_files); gc()

lfc_files <- list.files("Outputs/5_Predictability",
                        pattern = paste0("corLFCs_",MODE,".rds"),
                        full.names = T)

tissues <- c("Adipose-Subcutaneous","Artery-Tibial","Brain",
             "Esophagus-Mucosa","Lung","Muscle-Skeletal","Nerve-Tibial",
             "Skin-NotSunExposed(Suprapubic)","Skin-SunExposed(Lowerleg)",
             "Thyroid","WholeBlood")

for(tissue in tissues){
  
  lfc <- readRDS(lfc_files[grep(tissue,lfc_files,fixed = T)])
  genes <- c(names(head(sort(lfc[lfc > 1], decreasing = T),5)),
             names(head(sort(lfc[lfc < -1]),5)))
  
  young_expr <- readRDS(young_expr_files[grep(tissue,young_expr_files,
                                              fixed = T)])
  young_expr <- DataENSGToSymbol(young_expr)
  old_expr <- readRDS(old_expr_files[grep(tissue,old_expr_files,
                                          fixed = T)])
  old_expr <- DataENSGToSymbol(old_expr)
  young_pred <- readRDS(young_pred_files[grep(tissue,young_pred_files,
                                              fixed = T)])
  old_pred <- readRDS(old_pred_files[grep(tissue,old_pred_files,
                                          fixed = T)])
  
  pdf(paste0("Plots/5_Predictability/",tissue,"_tophits_",MODE,".pdf"),
      height = 7, width = 11)
  for(gene in genes){
    
    plot.data.o <- data.frame("Expression" = c(young_expr[gene,],old_expr[gene,]),
                              "Prediction" = c(young_pred[gene,],old_pred[gene,]),
                              "AgeGroup" = c(rep("Young",ncol(young_expr)),
                                             rep("Old",ncol(old_expr))))
    plot.data.o$AgeGroup <- factor(as.character(plot.data.o$AgeGroup),
                                   levels = c("Young","Old"))
    plot.data <- melt(plot.data.o)
    
    cors <- data.frame("AgeGroup" = c("Young","Old"),
                       "Cor" = c(readRDS(young_cor_files[grep(tissue,
                                                              young_cor_files,
                                                              fixed = T)])[gene],
                                 readRDS(old_cor_files[grep(tissue,
                                                            old_cor_files,
                                                            fixed = T)])[gene]))
    cors$Error <- log2(2-cors$Cor)
    cors$Cor <- round(cors$Cor,2)
    cors$Error <- round(cors$Error,2)
    
    plots <- list()
    
    plots[[1]] <- ggplot(plot.data) +
      geom_violin(aes(x = AgeGroup, y = value, fill = variable),
                  color = "black", draw_quantiles = 0.5, size = 1) +
      geom_point(aes(x = AgeGroup, y = value, fill = variable),
                 color = "black",
                 position = position_jitterdodge(dodge.width = .9)) +
      geom_text(data = cors, aes(x = AgeGroup, y = max(plot.data$value) + 0.5,
                                 label = paste0("error = ",Error))) +
      scale_fill_manual(values = c("Expression" = obs_color,
                                   "Prediction" = net_color),
                        labels = c("Expression" = "Observed",
                                   "Prediction" = "Predicted")) +
      xlab("Age") + ylab(paste0(gene," centered expr.")) +
      ggtitle(gsub("\\(","\n(",gsub("-"," - ",tissue)),
              sub = paste0("Error LFC = ",
                           as.character(round(lfc[gene],2)))) +
      labs(fill = NULL) + theme(text = element_text(size = 20),
                                legend.position = "bottom")
    
    plots[[2]] <- ggplot(plot.data.o) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_smooth(aes(x = Expression, y = Prediction, color = AgeGroup,
                      fill = AgeGroup), alpha = 0.2) +
      geom_point(aes(x = Expression, y = Prediction, color = AgeGroup)) +
      scale_color_manual(values = c("Young" = unname(age_palette["20-29"]),
                                    "Old" = unname(age_palette["50-59"]))) +
      scale_fill_manual(values = c("Young" = unname(age_palette["20-29"]),
                                   "Old" = unname(age_palette["50-59"]))) +
      xlab(paste0(gene," observed centered expr.")) +
      ylab(paste0(gene," predicted centered expr.")) +
      labs(color = "Age", fill = "Age") +
      theme(text = element_text(size = 20),
            legend.position = "bottom")
    
    print(ggarrange(plotlist = plots, nrow = 1, align = "hv"))
  }
  dev.off()
  
}

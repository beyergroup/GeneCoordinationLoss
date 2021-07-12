# Call poorly predicted genes

.libPaths <- c("/data/public/adesous1/GeneCorrelation/Resources/Rlibs/R-4.0.3", .libPaths())
library(reshape2)
library(ggplot2)

POORLY_PRED_COL = "tomato2"
net = "stabsel"
COR_FOLDER = paste0("Outputs/Human_Network/",net,"/Predictability/Tissue")

setwd("../")


# Load GTEx correlations

cor_files <- list.files(COR_FOLDER, pattern = "sampled_correlations.rds", full.names = T)
correlations <- sapply(cor_files, readRDS)
colnames(correlations) <- sapply(cor_files,
    function(x) strsplit(tail(strsplit(x, "/")[[1]], 1), "_")[[1]][1])
saveRDS(correlations, paste0("Outputs/Human_Network/",net,"/Predictability/Tissue/correlations_all.rds"))

# Call poorly predicted genes in cross tissue

poorly_predicted <- names(which(correlations[,"CrossTissue"] < COR_THRE)) # 4233 genes
poorly_predicted_all <- apply(correlations, 2, function(x) names(which(x < COR_THRE)))

saveRDS(poorly_predicted,
    paste0("Outputs/Human_Network/",net,"/Predictability/Tissue/poorly_predicted_crosstissue.rds"))
saveRDS(poorly_predicted_all,
    paste0("Outputs/Human_Network/",net,"/Predictability/Tissue/poorly_predicted_alltissues.rds"))


# Look at their correlations in individual tissues

plot.data <- melt(correlations)
plot.data$Group <- plot.data$Var1 %in% poorly_predicted
plot.data$Var2 <- gsub("(","\n(",plot.data$Var2,fixed=T)
plot.data$Var2 <- factor(as.character(plot.data$Var2),
    levels = c("CrossTissue","Adipose-Subcutaneous","Artery-Tibial","Brain",
        "Esophagus-Mucosa","Lung","Muscle-Skeletal","Nerve-Tibial",
        "Skin-NotSunExposed\n(Suprapubic)","Skin-SunExposed\n(Lowerleg)",
        "Thyroid","WholeBlood"))

pdf(paste0("Plots/Human_Network/",net,
      "/Predictability/poorly_predicted/CrossTissue_poorlypredicted_alltissues.pdf"),
    height = 8, width = 12)
ggplot(plot.data) +
  geom_violin(aes(y = value, x = Group, fill = Group), draw_quantiles = 0.5) +
  geom_hline(yintercept = 0.2, linetype = "dashed") +
  facet_wrap(~ Var2) + guides(fill = F) + xlab("") +
  ylab("Correlation coefficient") + ylim(c(-1,1)) +
  scale_fill_manual(values = c("TRUE" = POORLY_PRED_COL, "FALSE" = "darkgrey")) +
  scale_x_discrete(labels = c("TRUE" = "Poorly predicted", "FALSE" = "Others")) +
  theme(text = element_text(size = 20))
dev.off()


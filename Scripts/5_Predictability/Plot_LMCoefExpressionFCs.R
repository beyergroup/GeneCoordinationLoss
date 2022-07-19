# Plot delta LFCs vs expression LFCs

.libPaths("Resources/Rlibs/R-4.0.3/")
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
source("Scripts/functions.R")

# args = commandArgs(trailingOnly=TRUE)

predictability_files <- list.files("Outputs/5_Predictability",
                                   pattern = "_coefficients.rds",
                                   full.names = T)
expression_files <- list.files("Outputs/GTEx/AgeDE",
                               pattern = "_ageDE.rds", full.names = T)
# rmse_files <- list.files("Outputs/5_Predictability",
#                          pattern = "_rmse.rds", full.names = T)
cor_files <- list.files("Outputs/5_Predictability",
                        pattern = "_v7_sampled_centered_correlations.rds",
                        full.names = T)

well_predicted_genes <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")

tissues <- unique(sapply(predictability_files, function(x)
  strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1]))

plots <- list()
coef.plots <- list()
densities <- list()

for(tissue in tissues){
  
  ageDE <- readRDS(expression_files[grepl(tissue,expression_files,fixed =T)])
  ageDE <- limma::topTable(fit = ageDE, coef = "AgeGroupOld", number = nrow(ageDE))
  ageDE <- DataENSGToSymbol(ageDE, remove_dup = T)
  
  ageDP <- readRDS(predictability_files[grepl(tissue,predictability_files,fixed = T)])
  ageDP <- cbind(ageDP, "AgeFDR" = p.adjust(ageDP[,"AgePVal"]),
                 "ExprFDR" = p.adjust(ageDP[,"ExprPVal"]))
  
  # correct slope for the predictability of the gene across samples
  cor <- readRDS(cor_files[grepl(tissue,cor_files,fixed =T)])
  cor <- cor[rownames(ageDP)]
  ageDP <- cbind(ageDP, "CorrectedAgeSlope" = ageDP[,"AgeSlope"]*cor,
                 "CorrectedExprSlope" = ageDP[,"ExprSlope"]*cor)

  common <- intersect(rownames(ageDP),well_predicted_genes[[tissue]])

  plot.data <- data.frame("DE_LFC" = ageDP[common,"CorrectedExprSlope"],
                          "DE_pval" = ageDP[common,"ExprFDR"],
                          "DP_LFC" = ageDP[common,"CorrectedAgeSlope"],
                          "DP_pval" = ageDP[common,"AgeFDR"])

  coef.plots[[tissue]] <- ggplot() +
    geom_point(data = subset(plot.data, DE_pval > 0.05 & DP_pval > 0.1),
               aes(x = DE_LFC, y = DP_LFC), size = .5, color = "grey") +
    geom_point(data = subset(plot.data, DE_pval < 0.05),
               aes(x = DE_LFC, y = DP_LFC), color = "darkgrey") +
    geom_point(data = subset(plot.data, DP_pval < 0.1),
               aes(x = DE_LFC, y = DP_LFC), color = "orangered2") +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    ylim(c(-max(abs(plot.data$DP_LFC)),max(abs(plot.data$DP_LFC)))) +
    # ylim(c(-2,2)) +
    xlim(c(-max(abs(plot.data$DE_LFC)),max(abs(plot.data$DE_LFC)))) +
    xlab("Expression contribution") +
    ylab("Age contribution") +
    guides(alpha = F, color = F) + theme_classic() +
    ggtitle(gsub("(","\n(",tissue,fixed=T)) +
    theme(text = element_text(size = 20), title = element_text(size = 18))

  densities[[tissue]] <- ggplot(plot.data) +
    geom_density(aes(x = DP_LFC), size = .5, fill = "grey") +
    geom_vline(xintercept = 0) +
    xlim(c(-max(abs(plot.data$DP_LFC)),max(abs(plot.data$DP_LFC)))) +
    # xlim(c(-2,2)) +
    xlab("Error fold change") +
    guides(alpha = F, color = F) + theme_classic() +
    ggtitle(gsub("(","\n(",tissue,fixed=T)) +
    theme(text = element_text(size = 20), title = element_text(size = 18))

  
  common <- intersect(rownames(ageDP),rownames(ageDE))
  common <- intersect(common,well_predicted_genes[[tissue]])
  plot.data <- data.frame("DE_LFC" = ageDE[common,"logFC"],
                          "DE_pval" = ageDE[common,"adj.P.Val"],
                          "DP_LFC" = ageDP[common,"CorrectedAgeSlope"],
                          "DP_pval" = ageDP[common,"AgeFDR"])
  
  plots[[tissue]] <- ggplot() +
    geom_point(data = subset(plot.data, DE_pval > 0.05 & DP_pval > 0.1),
               aes(x = DE_LFC, y = DP_LFC), size = .5, color = "grey") +
    geom_point(data = subset(plot.data, DE_pval < 0.05),
               aes(x = DE_LFC, y = DP_LFC), color = "darkgrey") +
    geom_point(data = subset(plot.data, DP_pval < 0.1),
               aes(x = DE_LFC, y = DP_LFC), color = "orangered2") +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    ylim(c(-max(abs(plot.data$DP_LFC)),max(abs(plot.data$DP_LFC)))) +
    # ylim(c(-2,2)) +
    xlim(c(-max(abs(plot.data$DE_LFC)),max(abs(plot.data$DE_LFC)))) +
    xlab("Expression fold change") + ylab("Error fold change") +
    guides(alpha = F, color = F) + theme_classic() +
    ggtitle(gsub("(","\n(",tissue,fixed=T)) + 
    theme(text = element_text(size = 20), title = element_text(size = 18))
}

pdf(paste0("Plots/5_Predictability/predictability_GTEx_subsets_tissue_age_DE_nopoorlypred_LM_RMSE_pvalcutoff.pdf"),
    height = 10, width = 16)
print(ggarrange(plotlist = plots))
dev.off()

pdf(paste0("Plots/5_Predictability/predictability_GTEx_subsets_tissue_age_expr_nopoorlypred_LM_RMSE_pvalcutoff.pdf"),
    height = 10, width = 16)
print(ggarrange(plotlist = coef.plots))
dev.off()

pdf(paste0("Plots/5_Predictability/predictability_GTEx_subsets_tissue_age_nopoorlypred_LM_RMSE_pvalcutoff.pdf"),
    height = 8, width = 16)
print(ggarrange(plotlist = densities))
dev.off()


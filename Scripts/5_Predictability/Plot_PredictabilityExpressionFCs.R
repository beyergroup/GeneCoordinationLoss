# Plot delta LFCs vs expression LFCs

.libPaths("Resources/Rlibs/R-4.0.3/")
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
source("Scripts/functions.R")

# args = commandArgs(trailingOnly=TRUE)

predictability_files <- list.files(paste0("Outputs/Human_Network/stabsel/Predictability/AgeTissue"),
  pattern = "_sampled_ageDP.rds", full.names = T)
expression_files <- list.files("Outputs/GTEx/AgeDE",
  pattern = "_ageDE.rds", full.names = T)

well_predicted_genes <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")

tissues <- unique(sapply(predictability_files, function(x)
  strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1]))

plots <- list()
densities <- list()

for(tissue in tissues){
  
  ageDE <- readRDS(expression_files[grepl(tissue,expression_files,fixed =T)])
  ageDE <- limma::topTable(fit = ageDE, coef = "AgeGroupOld", number = nrow(ageDE))
  ageDE <- DataENSGToSymbol(ageDE, remove_dup = T)
  
  ageDP <- readRDS(predictability_files[grepl(tissue,predictability_files,fixed = T)])
  ageDP <- limma::topTable(fit = ageDP, coef = "AgeGroupOld", number = nrow(ageDP))
  
  common <- intersect(rownames(ageDP),rownames(ageDE))
  common <- intersect(common,well_predicted_genes[[tissue]])
  
  plot.data <- data.frame("DE_LFC" = ageDE[common,"logFC"], "DE_pval" = ageDE[common,"adj.P.Val"],
                          "DP_LFC" = ageDP[common,"logFC"], "DP_pval" = ageDP[common,"adj.P.Val"])
  
  plots[[tissue]] <- ggplot() +
    geom_point(data = subset(plot.data, DE_pval > 0.05 & DP_pval > 0.05),
               aes(x = DE_LFC, y = DP_LFC), size = .5, color = "grey") +
    geom_point(data = subset(plot.data, DE_pval < 0.05),
               aes(x = DE_LFC, y = DP_LFC), color = "darkgrey") +
    geom_point(data = subset(plot.data, DP_pval < 0.05),
               aes(x = DE_LFC, y = DP_LFC), color = "orangered2") +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    # ylim(c(-max(abs(plot.data$DP_LFC)),max(abs(plot.data$DP_LFC)))) +
    ylim(c(-2,2)) +
    xlim(c(-max(abs(plot.data$DE_LFC)),max(abs(plot.data$DE_LFC)))) +
    xlab("Expression fold change") + ylab("Error fold change") +
    guides(alpha = F, color = F) + theme_classic() +
    ggtitle(gsub("(","\n(",tissue,fixed=T)) + 
    theme(text = element_text(size = 20), title = element_text(size = 18))
  
  densities[[tissue]] <- ggplot(plot.data) +
    geom_density(aes(x = DP_LFC), size = .5, fill = "grey") +
    geom_vline(xintercept = 0) +
    # xlim(c(-max(abs(plot.data$DP_LFC)),max(abs(plot.data$DP_LFC)))) +
    xlim(c(-2,2)) +
    xlab("Error fold change") +
    guides(alpha = F, color = F) + theme_classic() +
    ggtitle(gsub("(","\n(",tissue,fixed=T)) + 
    theme(text = element_text(size = 20), title = element_text(size = 18))
}

pdf(paste0("Plots/5_Predictability/predictability_GTEx_subsets_tissue_age_DE_nopoorlypred.pdf"),
    height = 10, width = 16)
print(ggarrange(plotlist = plots))
dev.off()

pdf(paste0("Plots/5_Predictability/predictability_GTEx_subsets_tissue_age_nopoorlypred.pdf"),
    height = 8, width = 16)
print(ggarrange(plotlist = densities))
dev.off()

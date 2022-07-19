# Plot delta LFCs vs expression LFCs

.libPaths("Resources/Rlibs/R-4.0.3/")
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
source("Scripts/functions.R")

# args = commandArgs(trailingOnly=TRUE)

error_files <- list.files("Outputs/5_Predictability", pattern = "_errorLFCs.rds", full.names = T)
expression_files <- list.files("Outputs/GTEx/AgeDE", pattern = "_ageDE.rds", full.names = T)

well_predicted_genes <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")

tissues <- unique(sapply(error_files, function(x)
  strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1]))

# plots <- list()
# densities <- list()
plot.data <- data.frame()

for(tissue in tissues){
  
  ageDE <- readRDS(expression_files[grepl(tissue,expression_files,fixed =T)])
  ageDE <- limma::topTable(fit = ageDE, coef = "AgeGroupOld", number = nrow(ageDE))
  ageDE <- DataENSGToSymbol(ageDE, remove_dup = T)
  
  ageDP <- readRDS(error_files[grepl(tissue,error_files,fixed = T)])
  
  common <- intersect(names(ageDP),rownames(ageDE))
  common <- intersect(common,well_predicted_genes[[tissue]])
  
  plot.data <- rbind.data.frame(plot.data,
                                data.frame("DE_LFC" = ageDE[common,"logFC"],
                                           "DP_LFC" = ageDP[common],
                                           "Tissue" = tissue))
}  


pdf(paste0("Plots/5_Predictability/error_GTEx_subsets_tissue_age_DE_nopoorlypred.pdf"),
    height = 10, width = 16)
ggplot(plot.data) +
  geom_point(aes(x = DE_LFC, y = DP_LFC),
             size = .5, color = "grey") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ylim(c(-max(abs(plot.data$DP_LFC)),
         max(abs(plot.data$DP_LFC)))) +
  xlim(c(-max(abs(plot.data$DE_LFC)),
         max(abs(plot.data$DE_LFC)))) +
  xlab("Expression fold change") +
  ylab("Error fold change") +
  guides(alpha = F, color = F) +
  facet_wrap(~ Tissue) +
  theme_classic() +
  theme(text = element_text(size = 20),
        title = element_text(size = 18))
dev.off()

pdf(paste0("Plots/5_Predictability/error_GTEx_subsets_tissue_age_nopoorlypred.pdf"),
    height = 8, width = 16)
ggplot(plot.data) +
  geom_density(aes(x = DP_LFC), size = .5, fill = "grey") +
  geom_vline(xintercept = 0) +
  xlim(c(-max(abs(plot.data$DP_LFC)),max(abs(plot.data$DP_LFC)))) +
  xlab("Error fold change") +
  guides(alpha = F, color = F) +
  facet_wrap(~ Tissue) +
  theme_classic() +
  # ggtitle(gsub("(","\n(",tissue,fixed=T)) + 
  theme(text = element_text(size = 20), title = element_text(size = 18))
dev.off()

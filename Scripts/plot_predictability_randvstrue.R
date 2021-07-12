# Randomized vs true relationships

setwd("../")

library(reshape2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggplot2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggpubr, lib.loc = "Resources/Rlibs/R-4.0.3/")

type = "tissue"

correlations <- list()

for(net in c("stabsel","stabsel_randomized")){
  
  correlations[[net]] <- list()
  
  if(type == "age"){
    files <- list.files(paste0("Outputs/Human_Network/",net,"/Predictability/Age"),
                        pattern = "correlations.rds", full.names = T)
  } else if(type == "tissue"){
    files <- list.files(paste0("Outputs/Human_Network/",net,"/Predictability/Tissue"),
                        pattern = "correlations.rds", full.names = T)
  } else{
    files <- list.files(paste0("Outputs/Human_Network/",net,"/Predictability/AgeTissue"),
                        pattern = "correlations.rds", full.names = T)
  }
  
  for(file in files){
    
    prefix <- strsplit(tail(strsplit(file, "/")[[1]],1),"_sampled_correlations")[[1]][1]
    
    correlations[[net]][[prefix]] <- readRDS(file)
    
  }
}



# Plot

plot.data <- melt(correlations)
if(type = "agetissue"){
  plot.data$Tissue <- sapply(plot.data$L2, function(x) strsplit(x, "_")[[1]][1])
  plot.data$Tissue <- sapply(plot.data$Tissue, function(x) gsub("(","\n(",x, fixed = T))
  plot.data$Age <- sapply(plot.data$L2, function(x) strsplit(x, "_")[[1]][2])
} else if(type == "tissue"){
  plot.data$Tissue <- sapply(plot.data$L2, function(x) gsub("(","\n(",x, fixed = T))
}
rm(correlations); gc()

plots <- list()
if(type == "age"){
  
  for(net in unique(plot.data$L1)){
    plots[[net]] <- ggplot(subset(plot.data, L1 == net)) +
      geom_density(aes(x = value, color = L2), alpha = .4, size = 1) +
      scale_color_viridis_d(option = "magma", begin = .3, end = .7) +
      xlim(c(-1,1)) + geom_vline(xintercept = 0, linetype = "dashed") +
      xlab("Correlation coefficient") + ggtitle(net) +
      theme_classic() + theme(text = element_text(size = 20),
                              legend.position = "bottom", legend.title = element_blank(),
                              axis.title.y = element_blank(), title = element_text(size = 18))
  }
  
  pdf("Plots/Human_Network/Comparisons/predictability_GTEx_subsets_age.pdf")
  ggarrange(plotlist = plots, align = "hv", nrow = 2, common.legend = T, legend = "bottom")
  dev.off()
  
} else if(type == "tissue"){
  
  for(tissue in unique(plot.data$Tissue)[c(4,1:3,5:12)]){
    plots[[tissue]] <- eval(ggplot() +
      geom_density(data = subset(plot.data, (Tissue == tissue) & (L1 == "stabsel_randomized")),
        aes(x = value, color = L1, fill = L1)) +
      geom_density(data = subset(plot.data, (L1 == "stabsel") & (Tissue == tissue)),
        aes(x = value, color = L1, fill = L1), size = 3, alpha = 0.7) +
      scale_color_manual(values = c("stabsel" = "#1f73e0", "stabsel_randomized" = "darkgrey"),
          labels = c("stabsel" = "True relationships", "stabsel_randomized" = "Randomized relationships")) +
      scale_fill_manual(values = c("stabsel" = "#1f73e0", "stabsel_randomized" = "darkgrey"),
          labels = c("stabsel" = "True relationships", "stabsel_randomized" = "Randomized relationships")) +
      xlim(c(-1,1)) + geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
      xlab("Correlation coefficient") + ggtitle(tissue) +
      theme_classic() + theme(text = element_text(size = 20),
        legend.position = "bottom", legend.title = element_blank(),
        axis.title.y = element_blank(), title = element_text(size = 16)))
  }
  
  pdf("Plots/Human_Network/Comparisons/predictability_GTEx_subsets_truevsrand.pdf", height = 8, width = 13)
  ggarrange(plotlist = plots, align = "hv", common.legend = T, legend = "bottom")
  dev.off()
  
} else{
  
  for(net in unique(plot.data$Net)){
    plots[[net]] <- list()
    
    for(tissue in unique(plot.data$Tissue)){
      plots[[net]][[tissue]] <- ggplot(subset(plot.data, (Net == net) & (Tissue == tissue))) +
        geom_density(aes(x = value, color = Age), alpha = .4, size = 1) +
        scale_color_viridis_d(option = "magma", begin = .3, end = .7) +
        xlim(c(-1,1)) + geom_vline(xintercept = 0, linetype = "dashed") +
        xlab("Correlation coefficient") + ggtitle(tissue) +
        theme_classic() + theme(text = element_text(size = 20),
                                legend.position = "bottom", legend.title = element_blank(),
                                axis.title.y = element_blank(), title = element_text(size = 18))
    }
  }
  
  pdf("Plots/Human_Network/stabsel/Predictability/predictability_GTEx_subsets_tissue_age.pdf",
      height = 9, width = 16)
  ggarrange(plotlist = plots[[1]], align = "hv", common.legend = T, legend = "bottom")
  dev.off()
  
  pdf("Plots/Human_Network/stabsel_pcclasso_filter01/Predictability/predictability_GTEx_subsets_tissue_age.pdf",
      height = 9, width = 16)
  ggarrange(plotlist = plots[[2]], align = "hv", common.legend = T, legend = "bottom")
  dev.off()
}

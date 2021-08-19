# Plot correlations between network predictions and original data

args <- commandArgs(trailingOnly = T)
TYPE = args[1]
exclude_poorly_predicted = as.logical(args[2])

TYPE = "Tissue"
exclude_poorly_predicted = T

setwd("../")



correlations <- list()

for(net in c("stabsel","stabsel_pcclasso")){
  
  correlations[[net]] <- list()
  
  files <- list.files(paste0("Outputs/Human_Network/",net,"/Predictability/",TYPE),
                      pattern = "correlations.rds", full.names = T)
  
  for(file in files){
    
    prefix <- strsplit(tail(strsplit(file, "/")[[1]],1),"_sampled_correlations")[[1]][1]
    
    correlations[[net]][[prefix]] <- readRDS(file)
    
  }
}

if(exclude_poorly_predicted){
  
  for(net in c("stabsel","stabsel_pcclasso")){
    
    poorly_predicted <- readRDS(paste0("Outputs/Human_Network/",net,
      "/Predictability/Tissue/poorly_predicted_crosstissue.rds"))
    
    correlations[[net]] <- lapply(correlations[[net]],
      function(x) x[!(names(x) %in% poorly_predicted)])
  }
}

library(reshape2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggplot2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggpubr, lib.loc = "Resources/Rlibs/R-4.0.3/")

# Plot

plot.data <- reshape2::melt(correlations)
if(TYPE == "AgeTissue"){
  plot.data$Tissue <- sapply(plot.data$L2, function(x) strsplit(x, "_")[[1]][1])
  plot.data$Tissue <- sapply(plot.data$Tissue, function(x) gsub("(","\n(",x, fixed = T))
  plot.data$Age <- sapply(plot.data$L2, function(x) strsplit(x, "_")[[1]][2])
} else if(TYPE == "Tissue"){
  plot.data$Tissue <- sapply(plot.data$L2, function(x) gsub("(","\n(",x, fixed = T))
}
plot.data$Net <- "Full network"
plot.data$Net[grepl("pcclasso", plot.data$L1)] <- "No indirect effects"
rm(correlations); gc()

plots <- list()
if(TYPE == "Age"){
  
  for(net in unique(plot.data$Net)){
    plots[[net]] <- ggplot(subset(plot.data, Net == net)) +
      geom_density(aes(x = value, color = L2), alpha = .4, size = 1) +
      scale_color_viridis_d(option = "magma", begin = .3, end = .7) +
      xlim(c(-1,1)) + geom_vline(xintercept = 0, linetype = "dashed") +
      xlab("Correlation coefficient") + ggtitle(net) +
      theme_classic() + theme(text = element_text(size = 20),
        legend.position = "bottom", legend.title = element_blank(),
        axis.title.y = element_blank(), title = element_text(size = 18))
  }
  if(exclude_poorly_predicted){
    pdf("Plots/Human_Network/Comparisons/predictability_GTEx_subsets_age_nopoorlypred.pdf")
  } else{
    pdf("Plots/Human_Network/Comparisons/predictability_GTEx_subsets_age.pdf")
  }
  print(ggarrange(plotlist = plots, align = "hv", nrow = 2, common.legend = T, legend = "bottom"))
  dev.off()
  
} else if(TYPE == "Tissue"){
  
  for(tissue in unique(plot.data$Tissue)[c(4,1:3,5:12)]){
    plots[[tissue]] <- ggplot(subset(plot.data, Tissue == tissue)) +
      geom_density(aes(x = value, color = Net), alpha = .4, size = 1) +
      scale_color_viridis_d(option = "magma", begin = .3, end = .7) +
      xlim(c(-1,1)) + geom_vline(xintercept = 0, linetype = "dashed") +
      xlab("Correlation coefficient") + ggtitle(tissue) +
      theme_classic() + theme(text = element_text(size = 20),
                              legend.position = "bottom", legend.title = element_blank(),
                              axis.title.y = element_blank(), title = element_text(size = 18))
  }
  if(exclude_poorly_predicted){
    pdf("Plots/Human_Network/Comparisons/predictability_GTEx_subsets_nopoorlypred.pdf", height = 9, width = 16)
  } else{
    pdf("Plots/Human_Network/Comparisons/predictability_GTEx_subsets.pdf", height = 9, width = 16)
  }
  print(ggarrange(plotlist = plots, align = "hv", common.legend = T, legend = "bottom"))
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
  if(exclude_poorly_predicted){
    pdf("Plots/Human_Network/stabsel/Predictability/predictability_GTEx_subsets_tissue_age_nopoorlypred.pdf",
        height = 9, width = 16)
    print(ggarrange(plotlist = plots[[1]], align = "hv", common.legend = T, legend = "bottom"))
    dev.off()
    
    pdf("Plots/Human_Network/stabsel_pcclasso/Predictability/predictability_GTEx_subsets_tissue_age_nopoorlypredicted.pdf",
        height = 9, width = 16)
    print(ggarrange(plotlist = plots[[2]], align = "hv", common.legend = T, legend = "bottom"))
    dev.off()
  } else{
    pdf("Plots/Human_Network/stabsel/Predictability/predictability_GTEx_subsets_tissue_age.pdf",
        height = 9, width = 16)
    print(ggarrange(plotlist = plots[[1]], align = "hv", common.legend = T, legend = "bottom"))
    dev.off()
    
    pdf("Plots/Human_Network/stabsel_pcclasso/Predictability/predictability_GTEx_subsets_tissue_age.pdf",
        height = 9, width = 16)
    print(ggarrange(plotlist = plots[[2]], align = "hv", common.legend = T, legend = "bottom"))
    dev.off()
  }
  
}

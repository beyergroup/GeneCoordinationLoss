# Compute correlations between network predictions and original data

args <- commandArgs(trailingOnly = T)
PATTERN = args[1]
PRED_DIR = args[2]
PLOT_DIR = args[3]
MODE = args[4]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")
library(reshape2)
library(ggplot2)
library(ggpubr)


if(MODE == "spearman"){
  files <- list.files(PRED_DIR, pattern = paste0(PATTERN,"_correlations_",MODE))
  rand_files <- list.files(paste0(PRED_DIR,"/Randomized"),
                           pattern = paste0(PATTERN,"_correlations_",MODE))
} else{
  files <- list.files(PRED_DIR, pattern = paste0(PATTERN,"_correlations.rds"))
  rand_files <- list.files(paste0(PRED_DIR,"/Randomized"),
                           pattern = paste0(PATTERN,"_correlations.rds"))
}

plots <- list()

for(i in 1:length(files)){
  
  correlations <- ReadRDS(paste0(PRED_DIR,"/",files[i]))
  rand_correlations <- ReadRDS(paste0(PRED_DIR,"/Randomized/",rand_files[i]))
  
  t <- strsplit(files[i], split = "_")[[1]][1]
  t <- tissue_labels_dash_all[t]
  
  if(!all(names(correlations) == names(rand_correlations)))
    stop("Randomized correlation dimensions don't match non-randomized")
  
  plot.data <- data.frame("Original" = correlations,
                          "Random" = rand_correlations)
  plot.data <- melt(plot.data)
  
  plot.data$variable <- factor(as.character(plot.data$variable),
                               levels = c("Random","Original"))
  
  plots[[t]] <- ggplot() + 
    geom_density(data = subset(plot.data, variable == "Random"),
                 aes(x = value, color = variable, fill = variable)) +
    geom_density(data = subset(plot.data, variable == "Original"),
                 aes(x = value, color = variable, fill = variable), size = 3, alpha = 0.7) +
    xlab("Correlation coefficient") +
    scale_color_manual(values = c("Original" = "#ED2839", "Random" = "darkgrey"),
                       labels = c("Original" = "Blood network",
                                  "Random" = "Randomized blood network")) +
    scale_fill_manual(values = c("Original" = "#ED2839", "Random" = "darkgrey"),
                      labels = c("Original" = "Blood network",
                                 "Random" = "Randomized blood network")) +
    ggtitle(t) + xlim(c(-1,1)) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_classic() + theme(text = element_text(size = 20),
                            legend.title = element_blank(),
                            legend.position = "bottom",
                            axis.title.y = element_blank())
  
}


pdf(paste0(PLOT_DIR,"/performance_against_random.pdf"), height = 14, width = 18)
ggarrange(plotlist = plots[order(names(plots))], ncol = 5, nrow = 6, common.legend = T, legend = "bottom")
dev.off()

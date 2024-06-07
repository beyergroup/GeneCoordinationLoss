# Compute correlations between network predictions and original data

args <- commandArgs(trailingOnly = T)
BLOOD_DIR = args[1]
NET_DIR = args[2]
PLOT_DIR = args[3]
MODE = args[4]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")
library(reshape2)
library(ggplot2)
library(ggpubr)

PATTERN = "sampled_centered"

tissue_labels_dash <- c(tissue_labels_dash, "AllCells" = "All celltypes",
                        "BCells" = "B cells", "Monocytes" = "Monocytes",
                        "NKCells" = "NK cells", "TCells" = "T cells")

if(MODE == "spearman"){
  
  files <- list.files(NET_DIR, pattern = paste0(PATTERN,"_correlations_",MODE))
  rand_files <- list.files(paste0(NET_DIR,"/Randomized"),
                           pattern = paste0(PATTERN,"_correlations_",MODE))
  blood_files <- list.files(BLOOD_DIR, pattern = paste0(PATTERN,"_correlations_",MODE))
  rand_blood_files <- list.files(paste0(BLOOD_DIR,"/Randomized"),
                                 pattern = paste0(PATTERN,"_correlations_",MODE))
} else{
  
  files <- list.files(NET_DIR, pattern = paste0(PATTERN,"_correlations.rds"))
  rand_files <- list.files(paste0(NET_DIR,"/Randomized"),
                           pattern = paste0(PATTERN,"_correlations.rds"))
  blood_files <- list.files(BLOOD_DIR, pattern = paste0(PATTERN,"_correlations.rds"))
  rand_blood_files <- list.files(paste0(BLOOD_DIR,"/Randomized"),
                                 pattern = paste0(PATTERN,"_correlations.rds"))
}

plots <- list()

for(i in 1:length(files)){
  
  correlations <- ReadRDS(paste0(NET_DIR,"/",files[i]))
  rand_correlations <- ReadRDS(paste0(NET_DIR,"/Randomized/",rand_files[i]))
  blood_correlations <- ReadRDS(paste0(BLOOD_DIR,"/",blood_files[i]))
  rand_blood_correlations <- ReadRDS(paste0(BLOOD_DIR,"/Randomized/",rand_blood_files[i]))
  
  common <- intersect(names(correlations), names(blood_correlations))
  correlations <- correlations[common]
  rand_correlations <- rand_correlations[common]
  blood_correlations <- blood_correlations[common]
  rand_blood_correlations <- rand_blood_correlations[common]
  
  t <- strsplit(files[i], split = "_")[[1]][1]
  t <- tissue_labels_dash_all[t]
  
  if(!all(names(correlations) == names(rand_correlations)))
    stop("Randomized correlation dimensions don't match non-randomized")
  
  plot.data <- data.frame("Original" = correlations,
                          "Random" = rand_correlations,
                          "Blood" = blood_correlations,
                          "RandomBlood" = rand_blood_correlations)
  plot.data <- melt(plot.data)
  
  plot.data$variable <- factor(as.character(plot.data$variable),
                               levels = c("Random","RandomBlood",
                                          "Original","Blood"))
  plot.data$Tissue <- "Cross-tissue"
  plot.data$Tissue[grep("Blood",plot.data$variable)] <- "Blood"
  plot.data$Network <- "Original"
  plot.data$Network[grep("Rand",plot.data$variable)] <- "Random"
  
  plots[[t]] <- ggplot(plot.data) + 
    geom_density(aes(x = value, linetype = Network, color = Tissue,
                     size = Network)) +
    xlab("Correlation coefficient") +
    scale_color_manual(values = c("Blood" = "#ED2839", "Cross-tissue" =  "#1f73e0")) +
    scale_size_manual(values = c("Original" = 3, "Random" = 1)) +
    ggtitle(t) + xlim(c(-1,1)) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_classic() + theme(text = element_text(size = 20),
                            legend.title = element_blank(),
                            legend.position = "bottom",
                            axis.title.y = element_blank())
  
}


pdf(paste0(PLOT_DIR,"/performance_against_crosstissue_",MODE,".pdf"),
    height = 15, width = 22)
ggarrange(plotlist = plots[order(names(plots))], ncol = 5, nrow = 6,
          common.legend = T, legend = "bottom")
dev.off()

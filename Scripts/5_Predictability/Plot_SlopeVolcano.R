# Volcano plots for main tissues

SLOPE_DIR = "Outputs/5_Predictability/Age"
PLOT_DIR = "Plots/5_Predictability/Age"
COR_THRE = "well_predicted"

library(ggplot2)
source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")

slope_files <- list.files(SLOPE_DIR,
                          pattern = paste0("ageslope_",COR_THRE,".rds"),
                          full.names = T)
slopes <- sapply(slope_files, ReadRDS, simplify = F)
names(slopes) <- sapply(names(slopes),
                        function(x) strsplit(tail(strsplit(x, "/")[[1]],1),
                                             "_")[[1]][1])
rm(slope_files); gc()

# Restrict to tissues with signal
slopes <- slopes[names(slopes) %in% MAIN_TISSUES]

# Plot volcano

pdf(paste0(PLOT_DIR,"/ageslope_volcanoes.pdf"))
for(tissue in names(slopes)){
  
  plot.data <- as.data.frame(slopes[[tissue]])
  
  print(ggplot(plot.data) +
    geom_point(aes(x = Slope, y = -log10(pval),
                   color = -log10(pval))) +
    scale_colour_viridis_c(option = "inferno", begin = .1, end = .8) +
    ggtitle(tissue, "Regression volcano plot") +
    theme_classic() + theme(text = element_text(size = 20),
                            legend.position = "bottom"))
}
dev.off()

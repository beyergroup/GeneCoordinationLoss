# PCA of deltas (predicted - observed)

setwd("../")

library(ggplot2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(factoextra, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggpubr, lib.loc = "Resources/Rlibs/R-4.0.3/")

# Within-tissue patterns
files <- list.files("Outputs/Human_Network/stabsel/Predictability/AgeTissue",
  pattern = "sampled_delta.rds", full.names = T)

tissues <- sapply(files, function(x) strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1])

for(tissue in unique(tissues)){
  
  ff <- files[grepl(tissue,files,fixed=T)]
  data <- sapply(ff, readRDS, simplify = F)
  data <- do.call(cbind, data)
  
  plot.data <- data.frame("Age" = c(rep("20-29",31),rep("30-39",31),
    rep("40-49",31),rep("50-59",31),rep("60-69",31)))
  
  # PCA
  pca <- prcomp(t(data), center = F, scale. = F)
  plot.data <- cbind.data.frame(plot.data, pca$x)
  vars <- get_eigenvalue(pca)
  
  plots <- list()
  for(pc in 1:4){
    plots[[pc]] <- list()
    for(pcc in (pc+1):5){
       plots[[pc]][[pcc]] <- eval(ggplot(plot.data) +
        geom_point(aes(x = get(paste0("PC",pc)), y = get(paste0("PC",pcc)), color = Age)) +
        scale_color_viridis_d(option = "magma", begin = .3, end = .7) +
        # ggtitle(tissue) +
        theme_classic() +
        xlab(paste0("PC",pc," (",round(vars$variance.percent[pc],1),"%)")) +
        ylab(paste0("PC",pcc," (",round(vars$variance.percent[pcc],1),"%)")) +
        theme(text = element_text(size = 20), legend.position = "bottom",
              legend.title = element_blank(), title = element_text(size = 18)))
    }
    plots[[pc]] <- ggarrange(plotlist = plots[[pc]], nrow = 1, common.legend = T)
  }
  
  pdf(paste0("Plots/Human_Network/stabsel/Predictability/PCA_",tissue,"_deltas_age.pdf"),
    height = 14, width = 16)
  print(ggarrange(plotlist = plots, nrow = length(plots), common.legend = T, align = "hv"))
  dev.off()
}

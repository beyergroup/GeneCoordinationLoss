# Correlation between correlations and age

library(biomaRt, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(reshape2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggplot2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggpubr, lib.loc = "Resources/Rlibs/R-4.0.3/")

files <- list.files("Outputs/Human_Network/stabsel/Predictability/AgeTissue",
                    pattern = "sampled_correlations.rds", full.names = T)
tissues <- sapply(files, function(x) strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1])

ages <- c(25,35,45,55,65)

spearman_cors <- list()
for(tissue in unique(tissues)){
  
  cor <- sapply(files[grepl(tissue,files,fixed=T)], readRDS)
  sp <- apply(cor, 1, function(x){
    if(all(!is.na(x))){
      return(c(cor.test(x = x, y = ages, method = "spearman")$estimate,
               cor.test(x = x, y = ages, method = "spearman")$p.value))
    } else{
      return(NA)
    }})
  sp <- do.call(rbind,sp)
  rownames(sp) <- rownames(cor)
  colnames(sp) <- c("rho","p.value")
  
  spearman_cors[[tissue]] <- sp[order(abs(sp[,1]), decreasing = T),]
}
rm(sp,cor); gc()

saveRDS(spearman_cors, "Outputs/Human_Network/stabsel/Predictability/AgeTissue/agetrends.rds")


# p-values < 0.1
source('/data/public/adesous1/GeneCorrelation/Scripts/functions.R')
.libPaths(c("Resources/Rlibs/R-4.0.3",.libPaths()))
colors <- viridisLite::viridis(2, begin = .3, end = .7, option = "magma")

gobp_dims <- list("Adipose-Subcutaneous" = c(9,12), "Artery-Tibial" = c(5.2,10),
                  "Brain" = c(8,10), "Esophagus-Mucosa" = c(10,10), "Lung" = c(7,11),
                  "Muscle-Skeletal" = c(6.5,9.5), "Nerve-Tibial" = c(4,10),
                  "Skin-NotSunExposed(Suprapubic)" = c(5.5,9), "Skin-SunExposed(Lowerleg)" = c(4.2,11),
                  "Thyroid" = c(3,12), "WholeBlood" = c(7,11))

for(tissue in names(spearman_cors)){
  
  sp <- spearman_cors[[tissue]]
  sp <- sp[sp[,2] < 0.1,]
  sp <- sp[!is.na(sp[,1]),]
  
  # increased predictability
  up <- rownames(sp)[sp[,1] > 0]
  GObp_up <- GetGOEnrich(up, na.omit(rownames(sp)), "MF")
  up_plot <- PlotGOEnrich(GObp_up, colors[2], "Predictability increase with age")
  
  # decreased predictability
  down <- rownames(sp)[sp[,1] < 0]
  GObp_down <- GetGOEnrich(down, na.omit(rownames(sp)), "MF")
  down_plot <- PlotGOEnrich(GObp_down, colors[1], "Predictability decrease with age")
  
  # pdf(paste0("Plots/Human_Network/stabsel/Predictability/cor_trends_age/agetrends_",tissue,"_GOmf.pdf"),
  #   height = gobp_dims[[tissue]][1], width = gobp_dims[[tissue]][2])
  print(ggarrange(plotlist = list(up_plot,down_plot), nrow = 2, align = "hv",
            heights = c(nrow(GObp_up)+4,nrow(GObp_down)+4)))
  # dev.off()
  
  rm(sp,up,down,GObp_up,GObp_down,up_plot,down_plot); gc()
}

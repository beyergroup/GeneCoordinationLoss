# p-values < 0.1
source('/data/public/adesous1/GeneCorrelation/Scripts/functions.R')
.libPaths(c("Resources/Rlibs/R-4.0.3",.libPaths()))

colors <- viridisLite::viridis(2, begin = .3, end = .7, option = "magma")

gobp_dims <- list("Esophagus-Mucosa" = c(2,10), "Lung" = c(9,13),
                  "Skin-NotSunExposed(Suprapubic)" = c(2,10),
                  "Skin-SunExposed(Lowerleg)" = c(8,12), "WholeBlood" = c(2,10))

files <- list.files("Outputs/Human_Network/stabsel/Predictability/AgeTissue",
                    pattern = "_sampled_ageDP.rds", full.names = T)

for(file in files){
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  ageDP <- readRDS(file)
  ageDP <- limma::topTable(fit = ageDP, coef = "AgeGroupOld", number = nrow(ageDP))
  
  # increased predictability
  up <- rownames(subset(ageDP, (adj.P.Val < 0.05) & (logFC > 0)))
  GObp_up <- GetGOEnrich(up, rownames(ageDP), "BP")
  up_plot <- PlotGOEnrich(GObp_up, colors[2], "Increased deviation with age")
  
  # decreased predictability
  down <- rownames(subset(ageDP, (adj.P.Val < 0.05) & (logFC < 0)))
  GObp_down <- GetGOEnrich(down, rownames(ageDP), "BP")
  down_plot <- PlotGOEnrich(GObp_down, colors[1], "Decreased deviation with age")
  
  pdf(paste0("Plots/Human_Network/stabsel/Predictability/agetrends_",tissue,"_GObp.pdf"),
      height = gobp_dims[[tissue]][1], width = gobp_dims[[tissue]][2])
  up_plot
  # down_plot
  # ggarrange(plotlist = list(up_plot,down_plot), nrow = 2, align = "hv",
  #           heights = c(nrow(GObp_up)+4,nrow(GObp_down)+4))
  dev.off()
  
  rm(ageDP,up,down,GObp_up,GObp_down,up_plot,down_plot); gc()
}

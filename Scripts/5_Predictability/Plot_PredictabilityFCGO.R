# p-values < 0.1
source('/data/public/adesous1/GeneCorrelation/Scripts/functions.R')
.libPaths("Resources/Rlibs/R-4.0.3")


colors <- viridisLite::viridis(2, begin = .3, end = .7, option = "magma")

# gobp_dims <- list("Esophagus-Mucosa" = c(2,10), "Lung" = c(9,13),
#                   "Skin-NotSunExposed(Suprapubic)" = c(2,10),
#                   "Skin-SunExposed(Lowerleg)" = c(8,12), "WholeBlood" = c(2,10))
gobp_dims <- list("Brain" = c(10,10),
                  "Esophagus-Mucosa" = c(5,10), "Lung" = c(9,13),
                  "Skin-NotSunExposed(Suprapubic)" = c(2,10),
                  "Skin-SunExposed(Lowerleg)" = c(6,12), "WholeBlood" = c(2,10))

files <- list.files("Outputs/Human_Network/stabsel/Predictability/AgeTissue",
                    pattern = "_sampled_ageDP.rds", full.names = T)

well_predicted_genes <- readRDS("Outputs/5_Predictability/well_predicted_genes.rds")


for(file in files){
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  message("\n",tissue)
  
  ageDP <- readRDS(file)
  ageDP <- limma::topTable(fit = ageDP, coef = "AgeGroupOld", number = nrow(ageDP))
  ageDP <- ageDP[intersect(rownames(ageDP),well_predicted_genes[[tissue]]),]
  
  # increased predictability
  up <- rownames(subset(ageDP, (adj.P.Val < 0.05) & (logFC > 0)))
  GObp_up <- GetGOEnrich(up, rownames(ageDP), "BP")
  up_plot <- PlotGOEnrich(GObp_up, colors[2], "Increased deviation with age")
  
  # decreased predictability
  down <- rownames(subset(ageDP, (adj.P.Val < 0.05) & (logFC < 0)))
  GObp_down <- GetGOEnrich(down, rownames(ageDP), "BP")
  down_plot <- PlotGOEnrich(GObp_down, colors[1], "Decreased deviation with age")
  
  up_plot
  down_plot
  
  if(length(c(up,down)) > 0){
    message("Writing pdf")
    pdf(paste0("Plots/5_Predictability/predictability_FC_",tissue,"_GObp_nopoorlypred.pdf"),
        height = gobp_dims[[tissue]][1], width = gobp_dims[[tissue]][2])
    print(ggarrange(plotlist = list(up_plot,down_plot), nrow = 2, align = "hv",
                    heights = c(nrow(GObp_up)+4,nrow(GObp_down)+4)))
    dev.off()
  }
  rm(ageDP,up,down,GObp_up,GObp_down,up_plot,down_plot); gc()
}

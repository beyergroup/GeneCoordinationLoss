# Plot delta LFCs vs expression LFCs

.libPaths(c("Resources/Rlibs/R-4.0.3",.libPaths()))
library(topGO, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggplot2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggpubr, lib.loc = "Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")

args = commandArgs(trailingOnly=TRUE)
NET = args[1]
exclude_poorly_predicted = as.logical(args[2])


predictability_files <- list.files(paste0("Outputs/Human_Network/",NET,"/Predictability/AgeTissue"),
                                   pattern = "_sampled_ageDP.rds", full.names = T)
if(exclude_poorly_predicted)
  poorly_predicted <- readRDS(paste0("Outputs/Human_Network/",NET,
                                     "/Predictability/Tissue/poorly_predicted_crosstissue.rds"))

for(file in predictability_files){
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  ageDP <- readRDS(file)
  ageDP <- limma::topTable(fit = ageDP, coef = "AgeGroupOld", number = nrow(ageDP))
  
  if(exclude_poorly_predicted)
    ageDP <- ageDP[!(rownames(ageDP) %in% poorly_predicted),]
  
  if(sum((ageDP$adj.P.Val < 0.05) & (ageDP$logFC > 0)) >= 5){
    GO.bp.up <- GetGOEnrich(foreground = rownames(subset(ageDP, (adj.P.Val < 0.05) &
                                                          (logFC > 0))),
                            background = rownames(ageDP), go = "BP")
    GO.mf.up <- GetGOEnrich(foreground = rownames(subset(ageDP, (adj.P.Val < 0.05) &
                                                           (logFC > 0))),
                            background = rownames(ageDP), go = "MF")
    plots <- list()
    plots[[1]] <- PlotGOEnrich(GO.bp.up, col = "orangered2",
                               title = paste0("Increased deviation with age\n(",tissue,")"))
    plots[[2]] <- PlotGOEnrich(GO.mf.up, col = "orangered2",
                               title = paste0("Increased deviation with age\n(",tissue,")"))
    pdf(paste0("Plots/Human_Network/",NET,"/Predictability/predictability_GTEx_",tissue,"_GOup.pdf"),
        height = 4+(nrow(GO.bp.up)+nrow(GO.mf.up))/3, width = 11)
    print(ggarrange(plotlist = plots, heights = 4+c(nrow(GO.bp.up),nrow(GO.mf.up)),
                    ncol = 1, align = "v"))
    dev.off()
    rm(GO.bp.up,GO.mf.up,plots); gc()
  }
  if(sum((ageDP$adj.P.Val < 0.05) & (ageDP$logFC < 0)) >= 5){
    GO.bp.down <- GetGOEnrich(foreground = rownames(subset(ageDP, (adj.P.Val < 0.05) &
                                                             (logFC < 0))),
                              background = rownames(ageDP), go = "BP")
    GO.mf.down <- GetGOEnrich(foreground = rownames(subset(ageDP, (adj.P.Val < 0.05) &
                                                             (logFC < 0))),
                              background = rownames(ageDP), go = "MF")
    plots <- list()
    plots[[1]] <- PlotGOEnrich(GO.bp.down, col = "orangered2",
                               title = paste0("Decreased deviation with age\n(",tissue,")"))
    plots[[2]] <- PlotGOEnrich(GO.mf.down, col = "orangered2",
                               title = paste0("Decreased deviation with age\n(",tissue,")"))
    pdf(paste0("Plots/Human_Network/",NET,"/Predictability/predictability_GTEx_",tissue,"_GOdown.pdf"),
        height = 4+(nrow(GO.bp.down)+nrow(GO.mf.down))/3, width = 11)
    print(ggarrange(plotlist = plots, heights = 4+c(nrow(GO.bp.down),nrow(GO.mf.down)),
                    ncol = 1, align = "v"))
    dev.off()
    rm(GO.bp.down,GO.mf.down,plots); gc()
  }
  
}



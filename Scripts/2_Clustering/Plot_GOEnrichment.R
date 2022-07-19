# Plot GO enrichment in selected modules

args <- commandArgs(trailingOnly = T)
NET = args[1]
WDIR = args[2]
REPRESENTATION = args[3]
MODULE = as.character(args[4])

.libPaths("Resources/Rlibs/R-4.0.3/")
library(ggpubr)
source("Scripts/functions.R")
source(paste0("Scripts/",WDIR,"/params.R"))


load(paste0("Outputs/",WDIR,"/",NET,"_",REPRESENTATION,"_GO_weight01.RData"))

for(MODULE in names(GObp.list)[1:5]){
  plot.list <- list()
  plot.list[[1]] <- PlotGOEnrich(GObp.list[[MODULE]], col = COLORS[MODULE],
                                 title = paste0("GObp in module ",MODULE))
  plot.list[[2]] <- PlotGOEnrich(GOmf.list[[MODULE]], col = COLORS[MODULE],
                                 title = paste0("GOcc in module ",MODULE))
  plot.list[[3]] <- PlotGOEnrich(GOcc.list[[MODULE]], col = COLORS[MODULE],
                                 title = paste0("GOmf in module ",MODULE))
  
  heights <- c(nrow(GObp.list[[MODULE]]),
               nrow(GOmf.list[[MODULE]]),
               nrow(GOcc.list[[MODULE]])) + 3.5
  
  pdf(paste0("Plots/",WDIR,"/GO_weight01_module_",MODULE,".pdf"),
      height = PDF_DIMS[MODULE,"Height"],
      width = PDF_DIMS[MODULE,"Width"])
  ggarrange(plotlist = plot.list, ncol = 1, nrow = length(plot.list),
            heights = heights, align = "hv")
  dev.off()
  
}
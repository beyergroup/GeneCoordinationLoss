# Compute GO enrichment for top p-values

SLOPE_DIR = "Outputs/5_Predictability/Age/Age_MaxSubset"
PLOT_DIR = "Plots/5_Predictability/Age/Age_MaxSubset"
COR_THRE = "well_predicted"
INCL_70 = T

library(ggpubr)
source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")

if(INCL_70){
  slope_files <- list.files(SLOPE_DIR,
                            pattern = paste0("ageslope_",COR_THRE,".rds"),
                            full.names = T)
} else{
  slope_files <- list.files(SLOPE_DIR,
                            pattern = paste0("ageslope_",COR_THRE,"_no70.rds"),
                            full.names = T)
}
slopes <- sapply(slope_files, ReadRDS, simplify = F)
names(slopes) <- sapply(names(slopes),
                        function(x) strsplit(tail(strsplit(x, "/")[[1]],1),
                                             "_")[[1]][1])
rm(slope_files); gc()


# Select top p-values
GO <- list()
for(tissue in names(slopes)){
  
  message(tissue)
  
  data <- slopes[[tissue]][order(slopes[[tissue]][,"pval"]),]
  
  hits_up <- rownames(data)[1:N_TOP][data[1:N_TOP,"Slope"] > 0]
  message(length(hits_up)," positive slopes")
  hits_dw <- rownames(data)[1:N_TOP][data[1:N_TOP,"Slope"] < 0]
  message(length(hits_dw)," negative slopes")
  
  if(INCL_70){
    WriteRDS(list("Up" = hits_up, "Down" = hits_dw),
             paste0(SLOPE_DIR,"/",tissue,"_predictability_hits.rds"))
  } else{
    WriteRDS(list("Up" = hits_up, "Down" = hits_dw),
             paste0(SLOPE_DIR,"/",tissue,"_predictability_hits_no70.rds"))
  }

  GO[[tissue]] <- list("up" = list(), "dw" = list())
  GO[[tissue]][["up"]][["BP"]] <- GetGOEnrich(foreground = hits_up,
                                              background = rownames(data),
                                              go = "BP")
  GO[[tissue]][["up"]][["MF"]] <- GetGOEnrich(foreground = hits_up,
                                              background = rownames(data),
                                              go = "MF")
  GO[[tissue]][["dw"]][["BP"]] <- GetGOEnrich(foreground = hits_dw,
                                              background = rownames(data),
                                              go = "BP")
  GO[[tissue]][["dw"]][["MF"]] <- GetGOEnrich(foreground = hits_dw,
                                              background = rownames(data),
                                              go = "MF")
}

if(INCL_70){
  WriteRDS(GO, paste0(SLOPE_DIR,"/GOenrich_",N_TOP,".rds"))
} else{
  WriteRDS(GO, paste0(SLOPE_DIR,"/GOenrich_",N_TOP,"_no70.rds"))
}


# Plot GO enrichments
plots <- list("up" = list(), "dw" = list())

for(tissue in MAIN_TISSUES){
  
  if(nrow(GO[[tissue]][["up"]][["BP"]]) > 0){
    plots[["up"]][[tissue]] <- PlotGOEnrich(GO[[tissue]][["up"]][["BP"]][
      GO[[tissue]][["up"]][["BP"]][,"log2Enrichment"] > 2,],
                                            col = unname(age_palette["20-29"]),
                                            title = tissue_labels_dash[tissue])
  }
  if(nrow(GO[[tissue]][["dw"]][["BP"]]) > 0){
    plots[["dw"]][[tissue]] <- PlotGOEnrich(GO[[tissue]][["dw"]][["BP"]][
      GO[[tissue]][["dw"]][["BP"]][,"log2Enrichment"] > 2,],
                                            col = unname(age_palette["60-69"]),
                                            title = tissue_labels_dash[tissue])
  }
}

if(INCL_70){
  pdf(paste0(PLOT_DIR,"/GO_enrichments.pdf"), width = 20, height = 20)
  up <- ggarrange(ggarrange(plotlist = plots[['up']][c(1:2,6)], ncol = 1,
                      align = "v", heights = c(6,3,2.5)),
            ggarrange(plotlist = plots[['up']][3:5], ncol = 1,
                      align = "v", heights = c(9,7,7,5,1)),
            ncol = 2, widths = c(5,4))
  down <- ggarrange(ggarrange(plotlist = plots[['dw']][1:4], ncol = 1,
                      align = "v", heights = c(6,9,3.5,4)),
            ggarrange(plotlist = plots[['dw']][-c(1:4)], ncol = 1,
                      align = "v"),
            ncol = 2, widths = c(4,5))
  ggarrange(up, down, nrow = 2, labels = c("A","B"),
            font.label = list(size = 22))
} else{
  pdf(paste0(PLOT_DIR,"/GO_enrichments_no70.pdf"), width = 25, height = 15)
  ggarrange(plotlist = plots[['up']], nrow = 4, ncol = 2, align = "v")
  ggarrange(plotlist = plots[['dw']], nrow = 3, ncol = 2, align = "v")
}
dev.off()

# Plot heatmaps of module correlation, expression and predictability

# source("Scripts/functions.R")
library(reshape2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggplot2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(pheatmap, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(RColorBrewer, lib.loc = "Resources/Rlibs/R-4.0.3/")


TYPE = "expression"
# "variance"
# "predictability"
# "correlation"

NET = "stabsel"
# "stabsel_pcclasso"
METHOD = "greedy"
# "rwalk"

labels <- c("expression" = "Corrected expression", "variance" = "Expr. variance",
            "predictability" = "Predictability", "correlation" = "Within-module cor.")
method_labels <- c("greedy" = "Greedy", "rwalk" = "Random Walk")

for(TYPE in names(labels)){
  
  for(NET in c("stabsel","stabsel_pcclasso")){
    
    for(METHOD in c("greedy","rwalk")){
      
      m <- readRDS(list.files(paste0("Outputs/Human_Network/",NET,"/Topology/Modules"),
                              full.names = T, recursive = T, pattern = paste("mean",TYPE,METHOD, sep = "_")))
      plot.data <- melt(m, varnames = c("Module","Tissue"))
      plot.data <- subset(plot.data, !is.na(value))
      # 2 modules removed in stabsel
      
      # order by tissue-specificity of expression
      hc <- readRDS(paste0("Outputs/Human_Network/",NET,
                           "/Topology/Modules/TissueDE/module_clustering_bytissueDE_",METHOD,".rds"))
      plot.data$Module <- factor(as.character(plot.data$Module), levels = hc$labels[hc$order])
      clusters <- cutree(hc, h = max(hc$height/1.5))
      clusters <- clusters[hc$labels[hc$order]]
      
      # cluster annotation row
      palette <- rainbow(n = length(unique(clusters)))
      names(palette) <- unique(clusters)
      
      # absolute expr annotation row
      abs_expr <- readRDS(paste0("Outputs/Human_Network/",NET,
        "/Topology/Modules/TissueDE/CrossTissue_mean_absexpression_",METHOD,".rds"))
      abs_expr <- abs_expr[names(clusters)]
      
      scale_max <- switch(TYPE,
                          "expression" = max(abs(na.omit(subset(plot.data,
                                                                Tissue != "CrossTissue")$value))),
                          "variance" = max(abs(na.omit(subset(plot.data,
                                                              Tissue != "CrossTissue")$value))),
                          "predictability" = 1,
                          "correlation" = 1)
      
      # pdf(paste0("Plots/Human_Network/",NET,"/Topology/Modules/mean_",TYPE,"_modules_",METHOD,".pdf"),
      #     width = 20)
      # print(ggplot(subset(plot.data, Tissue != "CrossTissue")) +
      #   geom_tile(aes(y = Tissue, x= Module, fill = value)) +
      #   ggtitle("Mean module expression", subtitle = paste(METHOD,"algorithm")) +
      #   scale_fill_distiller(palette = "RdBu",
      #                        limits = c(-1,1)*scale_max) +
      #   labs(fill = labels[TYPE]) +
      #   theme(text = element_text(size = 20), legend.position = "bottom"))
      # dev.off()
      
      m <- m[,colnames(m) != "CrossTissue"]
      m <- t(m)
      
      m <- m[,colnames(m) %in% names(clusters)]
      clusters <- clusters[intersect(names(clusters), colnames(m))]
      m <- m[,hc$labels]
      
      pheatmap(m, cluster_rows = F, cluster_cols = hc, scale = "none",
               breaks = seq(from = -scale_max, to = scale_max, length.out = 501),
               color = colorRampPalette(c("blue","white","red"))(n = 501),
               cutree_cols = length(unique(clusters)), annotation_legend = F,
               annotation_col = data.frame("Expression" = abs_expr, "Group" = as.factor(clusters)),
               annotation_colors = list("Expression" = c("white","black"), "Group" = palette),
               annotation_names_col = F, cellwidth = 1.5, cellheight = 10,
               filename = paste0("Plots/Human_Network/",NET,"/Topology/Modules/mean_",
                                 TYPE,"_modules_",METHOD,"_heatmap.pdf"),
               show_colnames = F, main = paste0(labels[TYPE], " (",
                                                method_labels[METHOD]," modules)"))
    }
  }
}

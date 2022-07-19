# Figure 1

library(cowplot)

gtex.heat.plot <- readRDS("Outputs/0.5_Correlation/heatmap_CrossTissue_GTEx_ribosome_plot.rds")
TS.10X.heat.plot <- readRDS("Outputs/0.5_Correlation/heatmap_Merged_merged_10X_TS_ribosome_plot.rds")
TS.SS2.heat.plot <- readRDS("Outputs/0.5_Correlation/heatmap_Merged_merged_Smart-seq2_TS_ribosome_plot.rds")

gtex.lm.plot <- readRDS("Outputs/0.5_Correlation/trendlines_GTEx_RPLP0_RPLP1_ribosome_plot.rds")
TS.10X.lm.plot <- readRDS("Outputs/0.5_Correlation/trendlines_10X_TS_RPLP0_RPLP1_ribosome_plot.rds")
TS.SS2.lm.plot <- readRDS("Outputs/0.5_Correlation/trendlines_Smart-seq2_TS_RPLP0_RPLP1_ribosome_plot.rds")

pdf("Plots/Figures/Figure1.pdf", width = 10, height = 15)
plot_grid(gtex.heat.plot, gtex.lm.plot,
          TS.10X.heat.plot, TS.10X.lm.plot,
          TS.SS2.heat.plot, TS.SS2.lm.plot,
          nrow = 3, labels = "AUTO",
          align = "hv", axis = "bltr")
dev.off()

# to adjust:
# - make it square
# - match correlation scales between 3 heatmaps
# - match column and row heatmap order of TS data to GTEx ones
# - make heatmap column and row titles smaller
# - common legend
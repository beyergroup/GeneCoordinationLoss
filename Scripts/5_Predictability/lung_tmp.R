# Check top hits

library(ggplot2)
library(viridis)
source("Scripts/functions.R")

TISSUE = "Nerve-Tibial"

slopes <- ReadRDS(paste0("Outputs/5_Predictability/tmp/",
                         TISSUE,"_ageslope_0.75.rds"))
slopes <- as.data.frame(slopes)

# "Volcano" plot

pdf(paste0("Plots/5_Predictability/tmp/",TISSUE,"_volcano_pval.pdf"))
ggplot(slopes) +
  geom_point(aes(x = Slope, y = -log10(pval),
                 color = -log10(pval))) +
  scale_colour_viridis_c(option = "inferno", begin = .1, end = .8) +
  ggtitle(TISSUE, "Regression volcano plot") +
  theme_classic() + theme(text = element_text(size = 20),
                          legend.position = "bottom")
dev.off()



# Take top genes on each side

slopes$Score <- sign(slopes$Slope)*(-log10(slopes$pval))
slopes <- slopes[order(slopes$Score),]

N = 50

saveRDS(list("Down" = head(rownames(slopes),N),
             "Up" = tail(rownames(slopes),N)),
        paste0("Outputs/5_Predictability/tmp/",
               TISSUE,"_predictability_hits.rds"))


GO <- list("BP" = list(),
           "MF" = list())
GO$BP$down <- GetGOEnrich(head(rownames(slopes),N),
                          rownames(slopes),
                          go = "BP", algorithm = "weight01")
GO$MF$down <- GetGOEnrich(head(rownames(slopes),N),
                          rownames(slopes),
                          go = "MF", algorithm = "weight01")
GO$BP$up <- GetGOEnrich(tail(rownames(slopes),N),
                        rownames(slopes),
                        go = "BP", algorithm = "weight01")
GO$MF$up <- GetGOEnrich(tail(rownames(slopes),N),
                        rownames(slopes),
                        go = "MF", algorithm = "weight01")

PlotGOEnrich(GO$BP$down, age_palette["60-69"],
             "GObp enrichment in genes decreasing correlation with age")
ggsave(paste0(TISSUE,"_GObp_down_",N,".png"))

PlotGOEnrich(GO$BP$up, age_palette["20-29"],
             "GObp enrichment in genes increasing correlation with age")
ggsave(paste0(TISSUE,"_GObp_up_",N,".png"))

PlotGOEnrich(GO$MF$down, age_palette["60-69"],
             "GOmf enrichment in genes decreasing correlation with age")
ggsave(paste0(TISSUE,"_GOmf_down_",N,".png"))

PlotGOEnrich(GO$MF$up, age_palette["20-29"],
             "GOmf enrichment in genes increasing correlation with age")
ggsave(paste0(TISSUE,"_GOmf_up_",N,".png"))

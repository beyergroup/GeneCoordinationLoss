# Compare predictability of genes according to blood-specific network and cross-tissue cell line network

args <- commandArgs(trailingOnly = T)
BLOOD_DIR = "Outputs/tmp/NetTest/WellPredicted_TissueFilters" # args[1]
NET_DIR = "Outputs/5_Predictability/WellPredicted_TissueFilters" # args[2]
PLOT_DIR = "Plots/tmp/NetTest/WellPredicted_TissueFilters" # args[3]
MODE = "spearman" # args[4]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")
library(reshape2)
library(ggplot2)
library(ggpubr)

PATTERN = "sampled_centered"

# Load predictions from both networks

if(MODE == "spearman"){
  
  correlations <- ReadRDS(paste0(NET_DIR,"/WholeBlood_",PATTERN,
                                 "_correlations_",MODE,".rds"))
  rand_correlations <- ReadRDS(paste0(NET_DIR,"/Randomized/WholeBlood_",
                                      PATTERN,"_correlations_",MODE,".rds"))
  blood_correlations <- ReadRDS(paste0(BLOOD_DIR,"/WholeBlood_",PATTERN,
                                       "_correlations_",MODE,".rds"))
  rand_blood_correlations <- ReadRDS(paste0(BLOOD_DIR,"/Randomized/WholeBlood_",PATTERN,"_correlations_",MODE,".rds"))
  
} else{
  
  # not supported
}


# Which genes are not even predictable?

CrossTissue_unique <- names(correlations)[!(names(correlations) %in% names(blood_correlations))]
Blood_unique <- names(blood_correlations)[!(names(blood_correlations) %in% names(correlations))]

length(CrossTissue_unique)
# 12 592 genes (!) in the cross-tissue network but not in blood

length(Blood_unique)
# 109 (including some MT encoded)

# GO enrichment on Blood-unique genes
go.bp <- GetGOEnrich(Blood_unique, go = "BP",
                     union(names(correlations),names(blood_correlations)))
go.mf <- GetGOEnrich(Blood_unique, go = "MF",
                     union(names(correlations),names(blood_correlations)))
go.cc <- GetGOEnrich(Blood_unique, go = "CC",
                     union(names(correlations),names(blood_correlations)))

# Plot enrichments
pdf(paste0(PLOT_DIR,"/GOenrich_bloodnetunique.pdf"), height = 5, width = 9)
PlotGOEnrich(go.bp, col = "#ED2839", title = "Blood net only (GObp)")
PlotGOEnrich(go.mf, col = "#ED2839", title = "Blood net only (GOmf)")
PlotGOEnrich(go.cc, col = "#ED2839", title = "Blood net only (GOcc)")
dev.off()


# Performance differences

common <- intersect(names(correlations), names(blood_correlations))
# 4 032

plot.data <- data.frame("Blood" = blood_correlations[common],
                        "CrossTissue" = correlations[common],
                        # compute differences
                        "Delta" = blood_correlations[common] - correlations[common])

pdf(paste0(PLOT_DIR,"/networkperformance_scatterplot.pdf"),
    height = 5, width = 5)
ggplot(plot.data) +
  geom_point(aes(x = CrossTissue, y = Blood), size = 1, alpha = .7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0.2, linetype = "dashed", color = "#ED2839", size = 1) +
  geom_abline(slope = 1, intercept = -0.2, linetype = "dashed", color = "#1f73e0", size = 1) +
  ggtitle("Network performance on GTEx blood") +
  ylab("Blood network") + xlab("Cross tissue network") +
  theme_classic() + theme(text = element_text(size = 18))
dev.off()


blood_better <- rownames(plot.data)[plot.data$Delta > 0.2]
crosstissue_better <- rownames(plot.data)[plot.data$Delta < -0.2]

# GO enrichment on genes better predicted by the blood network
go.bp <- GetGOEnrich(blood_better, background = common, go = "BP")
go.mf <- GetGOEnrich(blood_better, background = common, go = "MF")
go.cc <- GetGOEnrich(blood_better, background = common, go = "CC")

# Plot enrichments
# png(paste0(PLOT_DIR,"/GObp_bloodnetbetter.png"), height = 700, width = 800)
pdf(paste0(PLOT_DIR,"/GObp_bloodnetbetter.pdf"),
    height = 10, width = 13)
PlotGOEnrich(go.bp, col = "#ED2839", title = "Blood net better (GObp)")
dev.off()
png(paste0(PLOT_DIR,"/GOmf_bloodnetbetter.png"), height = 200, width = 600)
PlotGOEnrich(go.mf, col = "#ED2839", title = "Blood net better (GOmf)")
dev.off()


# GO enrichment on genes better predicted by the cross tissue network
go.bp <- GetGOEnrich(crosstissue_better, background = common, go = "BP")
go.mf <- GetGOEnrich(crosstissue_better, background = common, go = "MF")
go.cc <- GetGOEnrich(crosstissue_better, background = common, go = "CC")

# Plot enrichments
# png(paste0(PLOT_DIR,"/GObp_crosstissuenetbetter.png"), height = 400, width = 850)
pdf(paste0(PLOT_DIR,"/GObp_crosstissuenetbetter.pdf"), height = 4, width = 13)
PlotGOEnrich(go.bp, col = "#1f73e0", title = "Cross-tissue net better (GObp)")
dev.off()
png(paste0(PLOT_DIR,"/GOmf_crosstissuenetbetter.png"), height = 230, width = 700)
PlotGOEnrich(go.mf, col = "#1f73e0", title = "Cross-tissue net better (GOmf)")
dev.off()
png(paste0(PLOT_DIR,"/GOcc_crosstissuenetbetter.png"), height = 200, width = 700)
PlotGOEnrich(go.cc, col = "#1f73e0", title = "Cross-tissue net better (GOcc)")
dev.off()


args = commandArgs(trailingOnly = T)
NET = "stabsel_filtered_trans_largestCC"
PATTERN = "Adjacency_undirected_weightssum_rownormalized"
EXPR_DIR = "Outputs/3_GTExDataPrep/Subset_Data/Max_Subset"
PRED_DIR = "Outputs/5_Predictability/Age/Age_MaxSubset"
CLSTR_DIR = "Outputs/2_Clustering"
OUTDIR = "Plots/5_Predictability/Age/Age_MaxSubset"


library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
library(gtools)
library(ggpubr)
library(ggplot2)
library(circlize)
source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")


# Network
net <- ReadRDS(paste0("Outputs/0_Preprocessing/",NET,"_network_Hs.rds"))

# Clusters
membership <- ReadRDS(paste0(CLSTR_DIR,"/",NET,"_",PATTERN,
                             "_reclustered_membership.rds"))

TISSUE = "WholeBlood"
GENE = "ALB"


hits <- ReadRDS(paste0(PRED_DIR,"/",TISSUE,"_predictability_hits.rds"))

expr_files <- list.files(EXPR_DIR, pattern = TISSUE)
expr_files <- expr_files[grep("sampled_centered",expr_files)]
data <- sapply(paste0(EXPR_DIR,"/",expr_files), ReadRDS)
names(data) <- sapply(expr_files, function(x) strsplit(x,"_")[[1]][2])
data <- lapply(data, DataENSGToSymbol, remove_dup = T)

predictors <- names(which(net[GENE,] != 0))
cor <- lapply(data,
              function(m) apply(m[intersect(rownames(m),predictors),],
                                1, cor, y = m[GENE,]))
cor <- do.call(cbind, cor)
h <- Heatmap(cor, cluster_columns = F, clustering_distance_rows = "pearson",
        column_title = "Age groups", column_title_side = "bottom",
        row_title = "Predictor genes",
        col = colorRamp2(breaks = c(-1,0,1),
                         colors =  c("#034780","white","#bf0a26")),
        heatmap_legend_param = list(title = paste("Correlation to",GENE),
                                    direction = "horizontal",
                                    title_position = "leftcenter",
                                    legend_width = unit(4,"cm"),
                                    title_gp = gpar(fontsize = 18),
                                    labels_gp = gpar(fontsize = 16)),
        row_names_gp = gpar(fontsize = 16), column_names_gp = gpar(fontsize = 16),
        column_title_gp = gpar(fontsize = 20))

h <- grid.grabExpr(draw(h, heatmap_legend_side = "bottom"))

sig <- apply(cor, 1, function(x) cor.test(x,c(25,35,45,55,65,75))$p.value)
genes <- names(which(sig < 0.05))

weighted_cor <- round(apply(cor, 2, mean, weights = net[GENE,rownames(cor)]),2)

plot.data <- melt(lapply(data, function(m) m[genes,,drop = F]))
plot.data$Target <- unlist(lapply(data, function(m) m[GENE,]))[plot.data$Var2]

g <- ggplot(plot.data) +
  geom_point(aes(y = value, x = Target), size = 1) +
  geom_smooth(aes(y = value, x = Target, color = Var1, group = Var1),
              alpha = .1, method = "lm") +
  facet_grid(Var1 ~ L1) +
  scale_color_discrete(guide = "none") +
  xlab(paste(GENE,"expr.")) + ylab("Neighbour expr.") +
  theme_classic() + theme(text = element_text(size = 20))


scat <- ReadRDS(paste0("Outputs/5_Predictability/Age/Age_MaxSubset/",
                       GENE,"_",TISSUE,"_scatterplot_predobs.rds"))

pdf(paste0(OUTDIR,"/",GENE,"_",TISSUE,".pdf"), width = 10, height = 10)
ggarrange(g, h, scat, NULL, nrow = 2, ncol = 2, widths = c(3,2),
          heights = c(3,2))
dev.off()

pdf(paste0("Plots/Figures/Parts/",GENE,"_",TISSUE,"_scatters.pdf"),
    height = 7,
    width = 10)
ggarrange(scat, nrow = 2, labels = "J", font.label = list(size = 22))
dev.off()

pdf(paste0("Plots/Figures/Parts/",GENE,"_",TISSUE,"_heatmap.pdf"),
    height = 3.5,
    # height = 5,
    width = 6)
ggarrange(h, labels = "I", font.label = list(size = 22))
dev.off()

# P-value distribution (main tissues)

args <- commandArgs(trailingOnly = TRUE)
SLOPE_DIR = args[1]
TISSUES = args[2]
INCL_70 = as.logical(args[3])



# SLOPE_DIR = "Outputs/5_Predictability/Age"
# PLOT_DIR = "Plots/5_Predictability/Age"
# COR_THRE = "well_predicted"
# TISSUES = "all"
# INCL_70 = T


library(reshape2)
library(ggplot2)
library(igraph)
library(cowplot)
library(ggpubr)
library(ComplexHeatmap)
source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")

if(INCL_70){
  bg_files <- list.files(SLOPE_DIR,
                         pattern = paste0("ageslopeBG_",COR_THRE,".rds"),
                         full.names = T)
} else{
  bg_files <- list.files(SLOPE_DIR,
                         pattern = paste0("ageslopeBG_",COR_THRE,"_no70.rds"),
                         full.names = T)
}
bg_slopes <- sapply(bg_files, ReadRDS)
names(bg_slopes) <- sapply(names(bg_slopes),
                           function(x) strsplit(tail(strsplit(x, "/")[[1]],1),
                                                "_")[[1]][1])
rm(bg_files); gc()

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

foreground.data <- melt(slopes)


# Slopes heatmap (panel C) ----------------------------------------------------

if(INCL_70){
  if(TISSUES == "all"){
    h <- ReadRDS(paste0(SLOPE_DIR,"/slope_heatmap_all.rds"))
  } else{
    h <- ReadRDS(paste0(SLOPE_DIR,"/slope_heatmap_main.rds"))
  }
} else{
  h <- ReadRDS(paste0(SLOPE_DIR,"/slope_heatmap_main_no70.rds"))
}
h <- grid.grabExpr(draw(h, heatmap_legend_side = "bottom"))


# Barplot of hits (Panel E) ---------------------------------------------------

if(!INCL_70){
  files <- list.files(SLOPE_DIR, pattern = "predictability_hits_no70.rds",
                      full.names = T)
} else{
  files <- list.files(SLOPE_DIR, pattern = "predictability_hits.rds",
                      full.names = T)
}

hits <- sapply(files, ReadRDS, simplify = F)
names(hits) <- sapply(files, function(x) strsplit(tail(strsplit(x,"/")[[1]],
                                                       1),"_")[[1]][1])
hit.data <- melt(hits)
if(TISSUES == "all"){
  hit.data$L1 <- factor(as.character(hit.data$L1),
                        levels = c("Thyroid",
                                   "Adipose-Visceral(Omentum)",
                                   "Esophagus-Mucosa",
                                   "Artery-Tibial",
                                   "WholeBlood",
                                   "Colon-Transverse",
                                   "Lung",
                                   "Muscle-Skeletal",
                                   "Brain",
                                   "Adipose-Subcutaneous",
                                   "Skin-SunExposed(Lowerleg)",
                                   "Skin-NotSunExposed(Suprapubic)",
                                   "Nerve-Tibial",
                                   "Breast-MammaryTissue",
                                   "Esophagus-Muscularis",
                                   "Testis"))
} else{
  hit.data <- subset(hit.data, L1 %in% MAIN_TISSUES)
  hit.data$L1 <- factor(as.character(hit.data$L1),
                        levels = c("Adipose-Visceral(Omentum)",
                                   "Artery-Tibial",
                                   "WholeBlood",
                                   "Brain",
                                   "Breast-MammaryTissue",
                                   "Esophagus-Mucosa",
                                   "Testis",
                                   "Thyroid"))
}

rm(files); gc()

b <- ggplot(hit.data) +
  geom_bar(aes(x = L1, fill = L2), position = "dodge",
           color = "white", width = 0.7) +
  coord_flip() + ylab("Gene number") + guides(fill = guide_legend(nrow = 2)) +
  scale_fill_manual(values =
                       c("Down" = unname(age_palette["60-69"]),
                         "Up" = unname(age_palette["20-29"])),
                     labels = c("Up" = "Increasing predictability w/ age",
                                "Down" = "Decreasing predictability w/ age"),
                     name = NULL) +
  scale_x_discrete(labels = tissue_labels_dash[as.character(unique(hit.data$L1))]) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() + theme(text = element_text(size = 20),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank(),
                          legend.position = "bottom")

if(INCL_70){
  WriteRDS(b, paste0(SLOPE_DIR,"/barplot_main.rds"))
} else{
  WriteRDS(b, paste0(SLOPE_DIR,"/barplot_main_no70.rds"))
}

pdf("Plots/5_Predictability/Age/Age_MaxSubset/hits_barplot.pdf", height = 4)
b
dev.off()


# P-value distribution plot (panel B) -----------------------------------------

plot.data <- melt(lapply(bg_slopes, function(m) m[,"pval",]),
                  value.name = "pval", varnames = c("Gene","Iteration"))

plot.ann <- data.frame("L1" = names(table(subset(plot.data,
                                                 Iteration == "1")$L1)),
                       "Count" = as.numeric(table(subset(plot.data,
                                                         Iteration == "1")$L1)))

# max pvalue at which tissue has hits
max_hit_pval <- c()
for(tissue in names(hits)){
  max_hit_pval <- c(max_hit_pval,
                    max(subset(foreground.data, (L1 == tissue) &
                                 (Var1 %in% unname(unlist(hits[[tissue]]))) &
                                 (Var2 == "pval"))$value))
}
max_hit_pval <- data.frame("Max" = max_hit_pval, "L1" = names(hits))

if(TISSUES == "all"){
  plot.data$L1 <- factor(as.character(plot.data$L1),
                         levels = unique(plot.data$L1)[c(1:3,16,4:15)])
  foreground.data$L1 <- factor(as.character(foreground.data$L1),
                               levels = unique(foreground.data$L1)[c(1:3,16,4:15)])
  max_hit_pval$L1 <- factor(as.character(max_hit_pval$L1),
                            levels = unique(max_hit_pval$L1)[c(1:3,16,4:15)])
  plot.ann$L1 <- factor(as.character(plot.ann$L1),
                        levels = unique(plot.ann$L1)[c(1:3,16,4:15)])
} else{
  plot.data <- subset(plot.data, L1 %in% MAIN_TISSUES)
  plot.data$L1 <- factor(as.character(plot.data$L1),
                         levels = MAIN_TISSUES)
  foreground.data <- subset(foreground.data, L1 %in% MAIN_TISSUES)
  foreground.data$L1 <- factor(as.character(foreground.data$L1),
                               levels = MAIN_TISSUES)
  max_hit_pval <- subset(max_hit_pval, L1 %in% MAIN_TISSUES)
  max_hit_pval$L1 <- factor(as.character(max_hit_pval$L1),
                            levels = MAIN_TISSUES)
  plot.ann <- subset(plot.ann, L1 %in% MAIN_TISSUES)
  plot.ann$L1 <- factor(as.character(plot.ann$L1),
                        levels = MAIN_TISSUES)
}



p <- ggplot(plot.data) +
  geom_density(aes(x = pval), color = NA, fill = "grey") +
  geom_density(data = subset(plot.data,
                             Iteration %in% sample(unique(plot.data$Iteration),
                                                   5)),
               aes(x = pval, group = Iteration, color = "Background"),
               size = 0.1) +
  geom_density(data = subset(foreground.data, Var2 == "pval"),
               aes(x = value, color = "Predictability ~ Age"), size = 1) +
  geom_vline(data = max_hit_pval,
             aes(xintercept = Max), linetype = "dashed") +
  geom_text(data = max_hit_pval,
            aes(x = 0.1, y = 0.15, label = round(Max,2))) +
  geom_text(data = subset(plot.ann, L1 %in% MAIN_TISSUES),
            aes(x = 0.5, y = 0.15, label = paste0(Count," genes"))) +
  xlab("Regression p-value") +
  scale_color_manual(values = c("Background" = "black",
                                "Predictability ~ Age" = "red")) +
  facet_wrap(~ L1,
             # nrow = 2,
             # ncol = 2,
             ncol = 3,
             labeller = labeller(L1 = tissue_labels_dash)) +
  theme_classic()  + theme(text = element_text(size = 20),
                           axis.text.x = element_text(angle = 40, hjust  = 1),
                           axis.title.y = element_blank(),
                           legend.title = element_blank(),
                           legend.position = "bottom")


pdf("Plots/5_Predictability/Age/Age_MaxSubset/pval_dists_well_predicted_main_hitpvals.pdf",
    width = 10)
p
dev.off()


# p.hist <- ggplot(plot.data) +
#   # geom_histogram(aes(x = pval, after_stat(ncount)), color = NA, fill = "grey", binwidth = 0.01,
#   #                boundary = 0) +
#   geom_histogram(data = subset(plot.data,
#                              Iteration %in% sample(unique(plot.data$Iteration),
#                                                    5)),
#                aes(x = pval, group = Iteration, color = "Background", after_stat(ncount)),
#                binwidth = 0.01, boundary = 0, alpha = 0.5, position = "identity",
#                fill = "transparent", linetype = "dashed") +
#   # geom_density(data = subset(plot.data,
#   #                            Iteration %in% sample(unique(plot.data$Iteration),
#   #                                                  5)),
#   #              aes(x = pval, group = Iteration, color = "Background"),
#   #              size = 0.1, boundary = 0) +
#   geom_histogram(data = subset(foreground.data, Var2 == "pval"),
#                aes(x = value, color = "Predictability ~ Age", after_stat(ncount)),
#                binwidth = 0.01, boundary = 0, alpha = 0.5, fill = "transparent",
#                linetype = "dashed", position = "identity") +
#   # geom_density(data = subset(foreground.data, Var2 == "pval"),
#   #              aes(x = value, color = "Predictability ~ Age"), size = 1) +
#   # geom_vline(data = max_hit_pval,
#   #            aes(xintercept = Max), linetype = "dashed") +
#   # geom_text(data = max_hit_pval,
#   #           aes(x = 0.1, y = 0.15, label = round(Max,2))) +
#   geom_text(data = subset(plot.ann, L1 %in% MAIN_TISSUES),
#             aes(x = 0.5, y = 0.15, label = paste0(Count," genes"))) +
#   xlab("Regression p-value") +
#   scale_color_manual(values = c("Background" = "black",
#                                 "Predictability ~ Age" = "red")) +
#   facet_wrap(~ L1,
#              # nrow = 2,
#              # ncol = 2,
#              ncol = 3,
#              labeller = labeller(L1 = tissue_labels_dash)) +
#   theme_classic()  + theme(text = element_text(size = 20),
#                            axis.text.x = element_text(angle = 40, hjust  = 1),
#                            axis.title.y = element_blank(),
#                            legend.title = element_blank(),
#                            legend.position = "bottom")
# 
# pdf("Plots/5_Predictability/Age/Age_MaxSubset/pval_dists_well_predicted_main_hist.pdf",
#     width = 10)
# p.hist
# dev.off()


p.slopes <- ggplot(plot.data) +
  geom_density(aes(x = Slope), color = NA, fill = "grey") +
  geom_density(data = subset(plot.data,
                             Iteration %in% sample(unique(plot.data$Iteration),
                                                   5)),
               aes(x = Slope, group = Iteration, color = "Background"),
               size = 0.1) +
  geom_density(data = subset(foreground.data, Var2 == "Slope"),
               aes(x = value, color = "Predictability ~ Age"), size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Regression slope") +
  scale_color_manual(values = c("Background" = "black",
                                "Predictability ~ Age" = "red")) +
  facet_wrap(~ L1,
             nrow = 2,
             # ncol = 2,
             labeller = labeller(L1 = tissue_labels_dash)) +
  theme_classic()  + theme(text = element_text(size = 20),
                           axis.text.x = element_text(angle = 40, hjust  = 1),
                           axis.title.y = element_blank(),
                           legend.title = element_blank(),
                           legend.position = "bottom")


if(INCL_70){
  WriteRDS(p, paste0(PLOT_DIR,"/pval_dist_density_main.rds"))
  # WriteRDS(p.slopes, paste0(PLOT_DIR,"/slope_dist_density_main.rds"))
} else{
  WriteRDS(p, paste0(PLOT_DIR,"/pval_dist_density_main_no70.rds"))
}

rm(plot.data,plot.ann,hits,hit.data,foreground.data,
   max_hit_pval,slopes,bg_slopes); gc()



# Scheme (panel A) ------------------------------------------------------------

set.seed(1)
dummy.data <- data.frame("Predicted" = jitter(1:10, amount = 0.8),
                         "Observed" = jitter(1:10, amount = 1))
pred.dummy <- ggplot(dummy.data) +
  geom_point(aes(x = Observed, y = Predicted), size = 3) +
  geom_smooth(aes(x = Observed, y = Predicted), color = "black",
              linetype = "dashed", method = "lm", se = F) +
  theme_classic() + theme(text = element_text(size = 16),
                          axis.text = element_blank(),
                          axis.ticks = element_blank())
set.seed(2)
dummy.data <- data.frame("Predicted" = jitter(1:10, amount = 5),
                         "Observed" = jitter(1:10, amount = 15))
unpred.dummy <- ggplot(dummy.data) +
  geom_point(aes(x = Observed, y = Predicted), size = 3) +
  geom_smooth(aes(x = Observed, y = Predicted), color = "black",
              linetype = "dashed", method = "lm", se = F) +
  theme_classic() + theme(text = element_text(size = 16),
                          axis.text = element_blank(),
                          axis.ticks = element_blank())

# graph <- make_graph(edges = c("A","B","A","C","D","A",
#                              "E","F","E","J","E","A"))
# V(graph)$color <- "black"
# pdf("Plots/Figures/Figure2_dummynet.pdf", width = 4, height = 4)
# plot(pred, vertex.label.dist = 3, size = 5,
#      vertex.label.family = "Helvetica",
#      edge.arrow.size = 1, edge.color = "black",
#      layout = layout_nicely)
# dev.off()

# Dummy scheme of predictability changes with age
dummy.data <- data.frame("Age" = c("20-29","30-39","40-49","50-59","60-69","70-79"),
                         "PredUp" = c(0.6,0.65,0.7,0.71,0.88,0.87),
                         "PredDown" = c(0.9,0.9,0.82,0.83,0.79,0.77),
                         "Nochange" = c(0.75,0.73,0.8,0.7,0.77,0.79))
dummy.data <- melt(dummy.data)
dummy.data$Age <- factor(as.character(dummy.data$Age),
                         levels = c("20-29","30-39","40-49","50-59","60-69","70-79"))
age.dummy <- ggplot(dummy.data) +
  geom_point(aes(x = Age, y = value, color = variable)) +
  geom_smooth(aes(x = Age, y = value, group = variable, color = variable),
              method = "lm", se = FALSE) +
  ylab("Predictability") +
  scale_color_manual(values = c("PredUp" = unname(age_palette["20-29"]),
                                "PredDown" = unname(age_palette["60-69"]),
                                "Nochange" = "grey"),
                     labels = c("PredUp" = "Increasing predictability w/ age",
                                "PredDown" = "Decreasing predictability w/ age",
                                "Nochange" = "No predictability change")) +
  guides(color = guide_legend(nrow = 3)) +
  theme_classic() + theme(text = element_text(size = 16),
                          axis.ticks.y = element_blank(),
                          axis.text.y = element_blank(),
                          legend.title = element_blank(),
                          legend.position = "bottom")

if(INCL_70){
  if(TISSUES == "all"){
    pdf("Plots/Figures/Parts/FigureS3_B.pdf", width = 25, height = 6)
    AB <- p.slopes
  } else{
    pdf("Plots/Figures/Parts/Figure2_AB.pdf", width = 17, height = 7)
    A <- ggarrange(NULL, ggarrange(pred.dummy,unpred.dummy,
                                   ncol = 2, nrow = 1, align = "hv"),
                   NULL, age.dummy, nrow = 4, heights = c(1,2,1,3))
    AB <- plot_grid(A,p, ncol = 2, rel_widths = c(5,14), labels = c("A","B"),
                    label_size = 22)
  }
} else{
  pdf("Plots/Figures/Parts/FigureS4_A.pdf", width = 6, height = 8)
  AB <- ggarrange(p, labels = "A", font.label = list(size = 22))
}
AB
dev.off()

if(TISSUES == "all"){
  pdf("Plots/Figures/Parts/FigureS3_A.pdf", width = 12, height = 5.7)
  p
  dev.off()
} else{
  if(INCLUDE_70){
    pdf("Plots/Figures/Parts/Figure2_B.pdf", width = 6, height = 8)
    p
    dev.off()
  }
}

# Slope heatmaps --------------------------------------------------------------

if(INCL_70){
  if(TISSUES == "all"){
    pdf("Plots/Figures/Parts/FigureS3_C.pdf", width = 18, height = 6)
    C <- plot_grid(h, labels = "C", label_size = 22)
  } else{
    pdf("Plots/Figures/Parts/Figure2_C.pdf", width = 8, height = 3)
    C <- plot_grid(h, labels = "C", label_size = 22)
  }
} else{
  pdf("Plots/Figures/Parts/FigureS4_B.pdf", width = 8, height = 3)
  C <- plot_grid(h, labels = "B", label_size = 22)
}
C
dev.off()


# Bar plots -------------------------------------------------------------------

if(INCL_70){
  if(TISSUES == "all"){
    pdf("Plots/Figures/Parts/FigureS3_D.pdf", width = 8, height = 6.5)
    D <- ggarrange(b, labels = "C", font.label = list(size = 22))
  } else{
    pdf("Plots/Figures/Parts/Figure2_E.pdf", width = 8, height = 4)
    D <- ggarrange(b, labels = "D", font.label = list(size = 22))
  }
} else{
  pdf("Plots/Figures/Parts/FigureS4_C.pdf", width = 8, height = 4)
  D <- ggarrange(b, labels = "C", font.label = list(size = 22))
}
D
dev.off()


# if(INCL_70){
#   if(TISSUES == "all"){
#     pdf("Plots/Figures/Parts/Figure2_DE_alltissues.pdf", width = 32, height = 6)
#     DE <- ggarrange(plotlist = list(b,scat), ncol = 2, widths = c(5,27),
#                     labels = c("D","E"), font.label = list(size = 22),
#                     common.legend = T, legend = "bottom", align = "h")
#   } else{
#     pdf("Plots/Figures/Parts/Figure2_DE.pdf", width = 19, height = 6)
#     DE <- ggarrange(plotlist = list(b,scat), ncol = 2, widths = c(5,14),
#                     labels = c("D","E"), font.label = list(size = 22),
#                     common.legend = T, legend = "bottom", align = "h")
#   }
# } else{
#   pdf("Plots/Figures/Parts/Figure2_DE_no70.pdf", width = 19, height = 6)
#   DE <- ggarrange(plotlist = list(b,scat), ncol = 2, widths = c(5,14),
#                   labels = c("D","E"), font.label = list(size = 22),
#                   common.legend = T, legend = "bottom", align = "h")
# }
# DE
# dev.off()


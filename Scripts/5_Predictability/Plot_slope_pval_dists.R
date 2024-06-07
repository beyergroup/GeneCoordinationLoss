

SLOPE_DIR = "Outputs/5_Predictability/Age"
PLOT_DIR = "Plots/5_Predictability/Age"
COR_THRE = "well_predicted"
TISSUES = "all"
INCL_70 = T

library(reshape2)
library(ggplot2)
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


# pvalue distribution
plot.data <- melt(lapply(bg_slopes, function(m) m[,"pval",]),
                  value.name = "pval", varnames = c("Gene","Iteration"))

plot.ann <- data.frame("L1" = names(table(subset(plot.data,
                                                 Iteration == "1")$L1)),
                       "Count" = as.numeric(table(subset(plot.data,
                                                         Iteration == "1")$L1)))
if(INCL_70){
  pdf(paste0(PLOT_DIR,"/pval_dists_",COR_THRE,".pdf"),
      width = 10, height = 8)
  # width = 6, height = 6)
} else{
  pdf(paste0(PLOT_DIR,"/pval_dists_",COR_THRE,"_no70.pdf"),
      width = 10, height = 8)
  # width = 6, height = 6)
}
ggplot(plot.data) +
  geom_density(aes(x = pval), color = NA, fill = "grey") +
  geom_density(data = subset(plot.data,
                             Iteration %in% sample(unique(plot.data$Iteration),
                                                   5)),
               aes(x = pval, group = Iteration, color = "Background"),
               size = 0.1) +
  geom_density(data = subset(foreground.data, Var2 == "pval"),
               aes(x = value, color = "Cor ~ Age"), size = 1) +
  geom_text(data = plot.ann,
            aes(x = 0.25, y = 0.15, label = paste0(Count," genes"))) +
  xlab("Regression p-value") +
  ggtitle("p-value distribution of Correlation ~ Age") +
  scale_color_manual(values = c("Background" = "black",
                                "Cor ~ Age" = "red")) +
  facet_wrap(~ L1) + theme(axis.title.y = element_blank(),
                           legend.title = element_blank(),
                           legend.position = "bottom")
dev.off()


p <- ggplot(subset(plot.data, L1 %in% MAIN_TISSUES)) +
  geom_density(aes(x = pval), color = NA, fill = "grey") +
  geom_density(data = subset(plot.data,
                             (Iteration %in% sample(unique(plot.data$Iteration),
                                                   5)) &
                               (L1 %in% MAIN_TISSUES)),
               aes(x = pval, group = Iteration, color = "Background"),
               size = 0.1) +
  geom_density(data = subset(foreground.data, (Var2 == "pval") &
                               (L1 %in% MAIN_TISSUES)),
               aes(x = value, color = "Predictability ~ Age"), size = 1) +
  geom_text(data = subset(plot.ann, L1 %in% MAIN_TISSUES),
            aes(x = 0.25, y = 0.15, label = paste0(Count," genes"))) +
  xlab("Regression p-value") +
  ggtitle("p-value distribution of Predictability ~ Age") +
  scale_color_manual(values = c("Background" = "black",
                                "Predictability ~ Age" = "red")) +
  facet_wrap(~ L1, nrow = 2, labeller = labeller(L1 = tissue_labels_dash)) +
  theme_classic()  + theme(text = element_text(size = 20),
                           axis.text.x = element_text(angle = 40, hjust  = 1),
                           axis.title.y = element_blank(),
                           legend.title = element_blank(),
                           legend.position = "bottom")
WriteRDS(p, "Outputs/5_Predictability/Age/Age_MaxSubset/pval_dist_density_main.rds")


rm(plot.data); gc()


# slope distribution

plot.data <- melt(lapply(bg_slopes, function(m) m[,"Slope",]),
                  value.name = "Slope", varnames = c("Gene","Iteration"))

if(TISSUES == "all"){
  plot.data$L1 <- factor(as.character(plot.data$L1),
                         levels = unique(plot.data$L1)[c(1:3,16,4:15)])
  foreground.data$L1 <- factor(as.character(foreground.data$L1),
                               levels = unique(foreground.data$L1)[c(1:3,16,4:15)])
} else{
  plot.data <- subset(plot.data, L1 %in% MAIN_TISSUES)
  plot.data$L1 <- factor(as.character(plot.data$L1),
                         levels = MAIN_TISSUES)
  foreground.data <- subset(foreground.data, L1 %in% MAIN_TISSUES)
  foreground.data$L1 <- factor(as.character(foreground.data$L1),
                               levels = MAIN_TISSUES)
}

if(INCL_70){
  if(TISSUES == "all"){
    pdf(paste0(PLOT_DIR,"/slope_dists_",COR_THRE,"_alltissues.pdf"),
        width = 25, height = 6)
  } else{
    pdf(paste0(PLOT_DIR,"/slope_dists_",COR_THRE,".pdf"),
        width = 10, height = 8)
    # width = 6, height = 6)
  }
} else{
  pdf(paste0(PLOT_DIR,"/slope_dists_",COR_THRE,"_no70.pdf"),
      width = 10, height = 8)
  # width = 6, height = 6)
}
ggplot(plot.data) +
  geom_density(aes(x = Slope), color = NA, fill = "grey") +
  geom_density(data = subset(plot.data,
                             Iteration %in% sample(unique(plot.data$Iteration),
                                                   5)),
               aes(x = Slope, group = Iteration, color = "Background"),
               size = 0.1) +
  geom_density(data = subset(foreground.data, Var2 == "Slope"),
               aes(x = value, color = "Predictability ~ Age"),
               linetype = "dashed",
               size = 1) +
  xlab("Regression slope") +
  # ggtitle("Slope distribution of Predictability ~ Age") +
  scale_color_manual(values = c("Background" = "black",
                                "Predictability ~ Age" = "red")) +
  facet_wrap(~ L1, nrow = 2, labeller = labeller(L1 = tissue_labels_dash)) +
  theme_classic()  + theme(text = element_text(size = 20),
                           axis.text.x = element_text(angle = 40, hjust  = 1),
                           axis.title.y = element_blank(),
                           legend.title = element_blank(),
                           legend.position = "bottom")
dev.off()

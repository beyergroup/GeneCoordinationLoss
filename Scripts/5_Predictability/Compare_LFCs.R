# Total LFCs vs cis and trans LFCs

source("Scripts/functions.R")
library(ggplot2)
library(ggpubr)

LFC_files <- list.files("Outputs/5_Predictability",
                        pattern = "corLFCs_YvsO_adj_spearman.rds")

tissues <- unique(sapply(LFC_files, function(x) strsplit(x, "_")[[1]][1]))

plots <- list("cis" = list(), "trans" = list(),
              "cistrans" = list())
densities <- list()

for(tissue in tissues){
  
  total_LFCs <- readRDS(paste0("Outputs/5_Predictability/",
                               tissue,"_corLFCs_YvsO_adj_spearman.rds"))
  cis_LFCs <- readRDS(paste0("Outputs/5_Predictability/",
                             tissue,"_cis_corLFCs_YvsO_adj_spearman.rds"))
  trans_LFCs <- readRDS(paste0("Outputs/5_Predictability/",
                               tissue,"_trans_corLFCs_YvsO_adj_spearman.rds"))
  
  plot.data <- data.frame("TotalLFC" = total_LFCs[names(cis_LFCs)],
                          "CisLFC" = cis_LFCs,
                          "TransLFC" = trans_LFCs[names(cis_LFCs)])
  
  plots[["cis"]][[tissue]] <- ggplot(plot.data) +
    geom_point(aes(x = CisLFC, y = TotalLFC),
               alpha = .3) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed") +
    geom_smooth(aes(x = CisLFC, y = TotalLFC),
                method = "lm") +
    ggtitle(gsub("(","\n(",tissue,fixed=T)) +
    theme_classic() + theme(text = element_text(size = 20))
  
  plots[["trans"]][[tissue]] <- ggplot(plot.data) +
    geom_point(aes(x = TransLFC, y = TotalLFC),
               alpha = .3) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed") +
    geom_smooth(aes(x = TransLFC, y = TotalLFC),
                method = "lm") +
    ggtitle(gsub("(","\n(",tissue,fixed=T)) +
    theme_classic() + theme(text = element_text(size = 20))
  
  plots[["cistrans"]][[tissue]] <- ggplot(plot.data) +
    geom_point(aes(x = CisLFC, y = TransLFC),
               alpha = .1) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed") +
    geom_smooth(aes(x = CisLFC, y = TransLFC),
                method = "lm") +
    ggtitle(gsub("(","\n(",tissue,fixed=T)) +
    theme_classic() + theme(text = element_text(size = 20))
  
  
  densities[[tissue]] <- ggplot(plot.data) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_density(aes(x = TotalLFC), color = net_color) +
    geom_density(aes(x = CisLFC), color = cis_color) +
    geom_density(aes(x = TransLFC), color = trans_color) +
    ggtitle(gsub("(","\n(",tissue,fixed=T)) +
    xlab("Error LFC") +
    theme_classic() + theme(text = element_text(size = 20),
                            axis.title.y = element_blank())
  
  rm(plot.data); gc()
}

pdf("Plots/5_Predictability/cis_trans_adj_spearman.pdf",
    width = 16, height = 10)
ggarrange(plotlist = densities, align = "hv")
ggarrange(plotlist = plots[["cis"]], align = "hv")
ggarrange(plotlist = plots[["trans"]], align = "hv")
ggarrange(plotlist = plots[["cistrans"]], align = "hv")
dev.off()


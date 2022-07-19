# Plot delta LFCs vs expression LFCs

.libPaths("Resources/Rlibs/R-4.0.3/")
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
source("Scripts/functions.R")

# MODE = "YvsO"
# MODE = "YvsO_adj"
MODE = "YvsO_adj_spearman"

# args = commandArgs(trailingOnly=TRUE)

error_files <- list.files("Outputs/5_Predictability",
                          pattern = paste0("_corLFCs_",MODE,".rds"),
                          full.names = T)
expression_files <- list.files("Outputs/GTEx/AgeDE",
                               pattern = "_ageDE.rds",
                               full.names = T)

well_predicted_genes <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")

tissues <- unique(sapply(error_files, function(x)
  strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1]))

# plots <- list()
densities <- list()
plot.data <- data.frame()

for(tissue in tissues){
  
  ageDE <- readRDS(expression_files[grepl(tissue,expression_files,fixed = T)])
  ageDE <- limma::topTable(fit = ageDE, coef = "AgeGroupOld", number = nrow(ageDE))
  ageDE <- DataENSGToSymbol(ageDE, remove_dup = T)
  
  ageDP <- readRDS(error_files[grepl(tissue,error_files,fixed = T)])
  
  common <- intersect(names(ageDP),rownames(ageDE))
  common <- intersect(common,well_predicted_genes[[tissue]])
  
  message(tissue)
  message("mean = ",round(mean(ageDP[common]),2))
  message("median = ",round(median(ageDP[common]),2),"\n")
  
  curr.plot.data <- data.frame("DE_LFC" = ageDE[common,"logFC"],
                               "DP_LFC" = ageDP[common],
                               "DE_pval" = ageDE[common,"adj.P.Val"],
                               "Tissue" = tissue)
  # curr.plot.data <- subset(curr.plot.data, abs(DP_LFC) < 5)
  plot.data <- rbind.data.frame(plot.data, curr.plot.data)
  
  d.data <- data.frame("x" = density(curr.plot.data$DP_LFC)$x,
                       "y" = density(curr.plot.data$DP_LFC)$y)
  d.data$Group <- "Center"
  d.data$Group[d.data$x <= -1] <- "LeftTail"
  d.data$Group[d.data$x >= 1] <- "RightTail"
  
  densities[[tissue]] <- ggplot() +
    geom_area(data = d.data,
              aes(x = x, y = y, group = Group, fill = Group),
              color = "black", outline.type = "full") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    annotate("text", x = mean(c(1,max(d.data$x))),
             y = max(d.data$y)/6, col = "#FB8500", size = 5,
             label = paste0("N = ",sum(curr.plot.data$DP_LFC >= 1))) +
    annotate("text", x = mean(c(-1,min(d.data$x))),
             y = max(d.data$y)/6, col = "#126782", size = 5,
             label = paste0("N = ",sum(curr.plot.data$DP_LFC <= -1))) +
    xlab("Error fold change") +
    scale_fill_manual(values = c("LeftTail" = "#126782",
                                 "Center" = "transparent",
                                 "RightTail" = "#FB8500")) +
    guides(fill = "none", color = "none") +
    theme_classic() +
    ggtitle(gsub("(","\n(",tissue,fixed=T)) + 
    theme(text = element_text(size = 20),
          title = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  
}


pdf(paste0("Plots/5_Predictability/corrYvsO_GTEx_subsets_tissue_age_DE_nopoorlypred_",
           MODE,".pdf"), height = 10, width = 16)
ggplot(plot.data) +
  # geom_point(aes(x = DE_LFC, y = DP_LFC),
  #            size = .5, color = "grey") +
  # geom_point(aes(x = DE_LFC, y = DP_LFC),
  #            color = "grey", alpha = .2) +
  geom_point(data = subset(plot.data, (DE_pval >= 0.05) &
                             (abs(DP_LFC) < 1)),
             aes(x = DE_LFC, y = DP_LFC), color = "grey", size = .5) +
  geom_point(data = subset(plot.data, DE_pval < 0.05),
             aes(x = DE_LFC, y = DP_LFC), color = "darkgrey") +
  geom_point(data = subset(plot.data, DP_LFC > 1),
             aes(x = DE_LFC, y = DP_LFC),
             color = age_palette["50-59"]) +
  geom_point(data = subset(plot.data, DP_LFC < -1),
             aes(x = DE_LFC, y = DP_LFC),
             color = age_palette["20-29"]) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  # ylim(c(-max(abs(plot.data$DP_LFC)),
  #        max(abs(plot.data$DP_LFC)))) +
  # coord_cartesian(ylim = c(-5,5)) +
  xlim(c(-max(abs(plot.data$DE_LFC)),
         max(abs(plot.data$DE_LFC)))) +
  xlab("Expression fold change") +
  ylab("Error fold change") +
  guides(alpha = F, color = F) +
  facet_wrap(~ Tissue) +
  theme_classic() +
  theme(text = element_text(size = 20),
        title = element_text(size = 18))
dev.off()

pdf(paste0("Plots/5_Predictability/corrYvsO_GTEx_subsets_tissue_age_nopoorlypred_",
           MODE,".pdf"), height = 8, width = 16)
ggplot(plot.data) +
  geom_density(aes(x = DP_LFC), size = .5, fill = "grey") +
  # geom_vline(xintercept = 0) +
  geom_vline(aes(xintercept = mean(DP_LFC))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(c(-max(abs(plot.data$DP_LFC)),max(abs(plot.data$DP_LFC)))) +
  xlab("Error fold change") +
  guides(alpha = F, color = F) +
  facet_wrap(~ Tissue) +
  theme_classic() +
  # coord_cartesian(xlim = c(-5,5)) +
  # ggtitle(gsub("(","\n(",tissue,fixed=T)) + 
  theme(text = element_text(size = 20), title = element_text(size = 18))
dev.off()


# remove outliers
ggplot(plot.data) +
  geom_density(aes(x = DP_LFC)) +
  facet_wrap(~ Tissue, scales = "free_x")

pdf(paste0("Plots/5_Predictability/corrYvsO_GTEx_subsets_tissue_age_Dvar_nopoorlypred_tails_",
           MODE,".pdf"), height = 8, width = 13)
ggarrange(plotlist = densities, align = "hv")
dev.off()


# Top error LFC per tissue
GObp.list <- list()
GOmf.list <- list()
GOcc.list <- list()

for(tissue in tissues){
  
  tmp <- subset(plot.data, Tissue == tissue)
  tmp <- tmp[order(tmp$DP_LFC),]
  
  # top decreasing errors
  GObp.list[[tissue]] <- list()
  GOmf.list[[tissue]] <- list()
  GOcc.list[[tissue]] <- list()
  GObp.list[[tissue]][["up"]] <- GetGOEnrich(foreground = head(rownames(tmp),100),
                                             background = rownames(tmp), go = "BP",
                                             algorithm = "weight01")
  GOmf.list[[tissue]][["up"]] <- GetGOEnrich(foreground = head(rownames(tmp),100),
                                             background = rownames(tmp), go = "MF",
                                             algorithm = "weight01")
  GOcc.list[[tissue]][["up"]] <- GetGOEnrich(foreground = head(rownames(tmp),100),
                                             background = rownames(tmp), go = "CC",
                                             algorithm = "weight01")
  
  # top increasing errors
  GObp.list[[tissue]][["down"]] <- GetGOEnrich(foreground = tail(rownames(tmp),100),
                                               background = rownames(tmp), go = "BP",
                                               algorithm = "weight01")
  GOmf.list[[tissue]][["down"]] <- GetGOEnrich(foreground = tail(rownames(tmp),100),
                                             background = rownames(tmp), go = "MF",
                                             algorithm = "weight01")
  GOcc.list[[tissue]][["down"]] <- GetGOEnrich(foreground = tail(rownames(tmp),100),
                                               background = rownames(tmp), go = "CC",
                                               algorithm = "weight01")
}
save(list = c("GObp.list","GOmf.list","GOcc.list"),
     file = paste0("Outputs/5_Predictability/corLFCs_top100_GOlists_",
                   MODE,".RData"))

for(tissue in tissues){
  
  plots <- list()
  if(nrow(GObp.list[[tissue]]$up) > 0){
    plots[[length(plots)+1]] <- PlotGOEnrich(GOenrich = GObp.list[[tissue]]$up,
                                             col = age_palette["50-59"],
                                             title = "Increased error with age (BP)")
  }
  if(nrow(GOmf.list[[tissue]]$up) > 0){
    plots[[length(plots)+1]] <- PlotGOEnrich(GOenrich = GOmf.list[[tissue]]$up,
                                             col = age_palette["50-59"],
                                             title = "Increased error with age (MF)")
  }
  if(nrow(GOcc.list[[tissue]]$up) > 0){
    plots[[length(plots)+1]] <- PlotGOEnrich(GOenrich = GOcc.list[[tissue]]$up,
                                             col = age_palette["50-59"],
                                             title = "Increased error with age (CC)")
  }
  if(length(plots) > 0){
    pdf(paste0("Plots/5_Predictability/",tissue,"_corLFCs_top100_up_GO_",MODE,
               ".pdf"),
        width = widths[[MODE]]["up",tissue],
        height = heights[[MODE]]["up",tissue])
    print(ggarrange(plotlist = plots, nrow = length(plots), align = "hv",
          heights = ratios[[MODE]][[tissue]][["up"]]))
    dev.off()
  }
  
  rm(plots); gc()
  
  plots <- list()
  if(nrow(GObp.list[[tissue]]$down) > 0){
    plots[[length(plots)+1]] <- PlotGOEnrich(GOenrich = GObp.list[[tissue]]$down,
                                             col = age_palette["20-29"],
                                             title = "Decreased error with age (BP)")
  }
  if(nrow(GOmf.list[[tissue]]$down) > 0){
    plots[[length(plots)+1]] <- PlotGOEnrich(GOenrich = GOmf.list[[tissue]]$down,
                                             col = age_palette["20-29"],
                                             title = "Decreased error with age (MF)")
  }
  if(nrow(GOcc.list[[tissue]]$down) > 0){
    plots[[length(plots)+1]] <- PlotGOEnrich(GOenrich = GOcc.list[[tissue]]$down,
                                             col = age_palette["20-29"],
                                             title = "Decreased error with age (CC)")
  }
  if(length(plots) > 0){
    pdf(paste0("Plots/5_Predictability/",tissue,"_corLFCs_top100_down_GO_",
               MODE,".pdf"),
        width = widths[[MODE]]["down",tissue],
        height = heights[[MODE]]["down",tissue])
    print(ggarrange(plotlist = plots, nrow = length(plots), align = "hv",
                    heights = ratios[[MODE]][[tissue]][["down"]]))
    dev.off()
  }
  rm(plots); gc()
  
}

heights <- list("YvsO_adj" = cbind("Adipose-Subcutaneous" = c("up" = 5,
                                                              "down" = 8),
                                   "Artery-Tibial" = c("up" = 3,
                                                       "down" = 4),
                                   "Brain" = c("up" = 2,
                                               "down" = NA),
                                   "WholeBlood" = c("up" = 2,
                                                    "down" = NA)),
                "YvsO_adj_spearman" = cbind("Adipose-Subcutaneous" = c("up" = 10,
                                                                       "down" = 8),
                                            "Artery-Tibial" = c("up" = 3,
                                                                "down" = 3),
                                            "Brain" = c("up" = 5,
                                                        "down" = 3),
                                            "Esophagus-Mucosa" = c("up" = NA,
                                                                   "down" = 3.5)))

widths <- list("YvsO_adj" = cbind("Adipose-Subcutaneous" = c("up" = 12,
                                                             "down" = 10),
                                  "Artery-Tibial" = c("up" = 10, "down" = 7),
                                  "Brain" = c("up" = 6, "down" = NA),
                                  "WholeBlood" = c("up" = 7, "down" = NA)),
               "YvsO_adj_spearman" = cbind("Adipose-Subcutaneous" = c("up" = 10,
                                                                      "down" = 12),
                                           "Artery-Tibial" = c("up" = 7,
                                                               "down" = 8),
                                           "Brain" = c("up" = 8,
                                                       "down" = 8.5),
                                           "Esophagus-Mucosa" = c("up" = NA,
                                                                  "down" = 8.5)))
ratios <- list("YvsO_adj" = list("Adipose-Subcutaneous" = list("up" = c(5,2),
                                                               "down" = c(5,2,2.5)),
                                 "Artery-Tibial" = list("up" = 1,
                                                        "down" = c(1,1)),
                                 "Brain" = list("up" = 1),
                                 "WholeBlood" = list("up" = 1)),
               "YvsO_adj_spearman" = list("Adipose-Subcutaneous" = list("up" = c(4,2,2),
                                                                        "down" = c(2.8,1)),
                                          "Artery-Tibial" = list("up" = 1,
                                                                 "down" = 1),
                                          "Brain" = list("up" = c(1,1),
                                                         "down" = 1),
                                          "Esophagus-Mucosa" = list("down" = 1)))


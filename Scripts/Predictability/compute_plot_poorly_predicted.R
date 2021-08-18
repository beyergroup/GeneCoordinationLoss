# Call poorly predicted genes

.libPaths(c("/data/public/adesous1/GeneCorrelation/Resources/Rlibs/R-4.0.3", .libPaths()))
library(topGO, lib.loc = "/data/public/adesous1/GeneCorrelation/Resources/Rlibs/R-4.0.3")
library(reshape2)
library(ggplot2)
source("functions.R")

POORLY_PRED_COL = "tomato2"
POORLY_PRED_H_COL = "tomato3"
POORLY_PRED_L_COL = "lightcoral"
COR_THRE = 0.2
MEAN_THRE = 8
LOW_MEAN_THRE = 2

net = "stabsel"
NET_FILE = "Data/Networks/Human/stabsel_network_Hs_filtered.rds"
# net = "stabsel_pcclasso_filter01"
# NET_FILE = "Data/Networks/Human/stabsel_pcclasso_network_Hs_filtered.rds"

COR_FOLDER = paste0("Outputs/Human_Network/",net,"/Predictability/Tissue")
palette <- c("Poorly Predicted\n(High expr.)" = POORLY_PRED_H_COL,
             "Poorly Predicted\n(Low expr.)" = POORLY_PRED_L_COL,
             "Others" = "darkgrey")

setwd("../")


# Load GTEx correlations

cor_files <- list.files(COR_FOLDER, pattern = "sampled_correlations.rds", full.names = T)
correlations <- sapply(cor_files, readRDS)
colnames(correlations) <- sapply(cor_files,
    function(x) strsplit(tail(strsplit(x, "/")[[1]], 1), "_")[[1]][1])
saveRDS(correlations, paste0("Outputs/Human_Network/",net,"/Predictability/Tissue/correlations_all.rds"))


# Call poorly predicted genes in cross tissue

poorly_predicted <- names(which(correlations[,"CrossTissue"] < COR_THRE)) # 4217 genes
poorly_predicted_all <- apply(correlations, 2, function(x) names(which(x < COR_THRE)))

saveRDS(poorly_predicted,
    paste0("Outputs/Human_Network/",net,"/Predictability/Tissue/poorly_predicted_crosstissue.rds"))
saveRDS(poorly_predicted_all,
    paste0("Outputs/Human_Network/",net,"/Predictability/Tissue/poorly_predicted_alltissues.rds"))


# Load GTEx data
library(data.table)
GTEx_real <- fread("/data/public/adesous1/scDropImp/networkinference/analyses/paper/NetworkQuality_analysis/GTEx_DESeq2/GTEx_DESeq2_GEinput.txt",
                   header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)
GTEx_real <- data.frame(GTEx_real,row.names = 1)
GTEx_real <- as.matrix(GTEx_real[,-(1:2)]) # dim 35488 17382
GTEx_real <- log2(GTEx_real + 1)
# Reduce to common targets
GTEx_real <- GTEx_real[intersect(rownames(correlations),rownames(GTEx_real)),]
correlations <- correlations[rownames(GTEx_real),]
poorly_predicted <- poorly_predicted[poorly_predicted %in% rownames(GTEx_real)]

info <- data.frame("GTExMeans" = rowMeans(GTEx_real),
                   "GTExVars" = apply(GTEx_real,1,var))


# Expression across tissues
tissues <- readRDS("/data/public/adesous1/scDropImp/networkinference/analyses/paper/NetworkQuality_analysis/tissues.rds")
tissues <- droplevels(tissues)
tissue_means <- sapply(levels(tissues), function(t) rowMeans(GTEx_real[,as.character(tissues) == t]))
tissue_expressed <- rowSums(tissue_means > LOW_MEAN_THRE)
info$TissuesExpressed <- tissue_expressed[rownames(info)]

rm(GTEx_real,tissues,tissue_means,tissue_expressed); gc()


# Load network
network_matrix <- as.matrix(readRDS(NET_FILE))
network_matrix <- network_matrix[rownames(info),]

info$OriginalPredNumber <- rowSums(network_matrix[,] != 0)[rownames(info)]
network_matrix <- network_matrix[,intersect(colnames(network_matrix),rownames(info))]
info$QuantifiedPredNumber <- rowSums(network_matrix[,] != 0)[rownames(info)]
info$LowlyExpressedPredNumber <- sapply(rownames(info), function(g){
  predictors <- names(which(network_matrix[g,] != 0))
  return(sum(info$GTExMeans[rownames(info) %in% predictors] < LOW_MEAN_THRE))
})[rownames(info)]
info$WellExpressedPredNumber <- sapply(rownames(info), function(g){
  predictors <- names(which(network_matrix[g,] != 0))
  return(sum(info$GTExMeans[rownames(info) %in% predictors] >= LOW_MEAN_THRE))
})[rownames(info)]
info$WellExpressedPredFraction <- info$WellExpressedPredNumber/info$OriginalPredNumber
info$LowlyExpressedPredFraction <- info$LowlyExpressedPredNumber/info$OriginalPredNumber
rm(network_matrix); gc()


# Split poorly predicted into highly and lowly expressed
poorly_predicted_H <- poorly_predicted[info[poorly_predicted,"GTExMeans"] > MEAN_THRE]
poorly_predicted_L <- poorly_predicted[info[poorly_predicted,"GTExMeans"] <= MEAN_THRE]
saveRDS(poorly_predicted_H,
        paste0("Outputs/Human_Network/",net,"/Predictability/Tissue/poorly_predicted_crosstissue_H.rds"))
saveRDS(poorly_predicted_L,
        paste0("Outputs/Human_Network/",net,"/Predictability/Tissue/poorly_predicted_crosstissue_L.rds"))
info <- info[!is.na(correlations[,"CrossTissue"][rownames(info)]),]
info$Group <- "Others"
info$Group[rownames(info) %in% poorly_predicted_H] <- "Poorly Predicted\n(High expr.)"
info$Group[rownames(info) %in% poorly_predicted_L] <- "Poorly Predicted\n(Low expr.)"
# info$Group <- factor(as.character(info$Group), levels = c("Poorly Predicted\n(High expr.)",
#                                                           "Poorly Predicted\n(Low expr.)",
#                                                           "Others"))


# Look at their correlations in individual tissues

plot.data <- melt(correlations)
plot.data$Group <- "Others"
plot.data$Group[plot.data$Var1 %in% poorly_predicted_H] <- "Poorly Predicted\n(High expr.)"
plot.data$Group[plot.data$Var1 %in% poorly_predicted_L] <- "Poorly Predicted\n(Low expr.)"
plot.data$Group <- factor(as.character(plot.data$Group), levels = c("Poorly Predicted\n(High expr.)",
                                                                    "Poorly Predicted\n(Low expr.)",
                                                                    "Others"))
plot.data$Var2 <- gsub("(","\n(",plot.data$Var2,fixed=T)
plot.data$Var2 <- factor(as.character(plot.data$Var2),
    levels = c("CrossTissue","Adipose-Subcutaneous","Artery-Tibial","Brain",
        "Esophagus-Mucosa","Lung","Muscle-Skeletal","Nerve-Tibial",
        "Skin-NotSunExposed\n(Suprapubic)","Skin-SunExposed\n(Lowerleg)",
        "Thyroid","WholeBlood"))

pdf(paste0("Plots/Human_Network/",net,
      "/Predictability/poorly_predicted/CrossTissue_poorlypredicted_alltissues_H_L.pdf"),
    height = 8, width = 10.5)
ggplot(plot.data) +
  geom_violin(aes(y = value, x = Group, fill = Group), draw_quantiles = 0.5) +
  geom_hline(yintercept = 0.2, linetype = "dashed") +
  facet_wrap(~ Var2) + guides(fill = F) + xlab("") +
  ylab("Correlation coefficient") + ylim(c(-1,1)) +
  scale_fill_manual(values = palette) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
dev.off()

rm(plot.data); gc()


# Look at features of poorly and well predicted genes

plots <- list()

plots$Mean <- ggplot() +
  geom_density(data = subset(info, Group == "Others"),
               aes(x = GTExMeans), fill = palette["Others"],
               color = palette["Others"]) +
  geom_density(data = subset(info, Group != "Others"),
    aes(x = GTExMeans, fill = Group,
        color = Group), alpha = 0.6) +
  xlab("Mean expression across GTEx samples") +
  guides(fill = F, color = F, alpha = F) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) +
  theme_classic() + theme(text = element_text(size = 20),
                          axis.ticks.y = element_line(color = "transparent"),
                          axis.text.y = element_text(color = "transparent"),
                          axis.title.y = element_blank())

plots$Var <- ggplot() +
  geom_density(data = subset(info, Group == "Others"),
               aes(x = GTExVars), fill = palette["Others"],
               color = palette["Others"]) +
  geom_density(data = subset(info, Group != "Others"),
               aes(x = GTExVars, fill = Group,
                   color = Group), alpha = 0.6) +
  xlab("Variance across GTEx samples") +
  guides(fill = F, color = F, alpha = F) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) +
  theme_classic() + theme(text = element_text(size = 20),
                          axis.ticks.y = element_line(color = "transparent"),
                          axis.text.y = element_text(color = "transparent"),
                          axis.title.y = element_blank())

info$Group <- factor(as.character(info$Group),
                     levels = c("Others",
                                "Poorly Predicted\n(Low expr.)",
                                "Poorly Predicted\n(High expr.)"))

plots$PredNumber <- ggplot(info) +
  geom_bar(aes(x = OriginalPredNumber, fill = Group)) +
  scale_fill_manual(values = palette) +
  xlab("Number of predictors") + ylab("# genes") +
  theme_classic() + theme(text = element_text(size = 20),
                          legend.title = element_blank())

plots$Tissues <- ggplot(info) +
  geom_bar(aes(x = TissuesExpressed, fill = Group)) +
  scale_fill_manual(values = palette) +
  xlab("# tissues where gene is expressed") + ylab("# genes") +
  coord_cartesian(ylim = c(0,3000)) +
  theme_classic() + theme(text = element_text(size = 20),
                          legend.title = element_blank())

plots$Frac <- ggplot(info) +
  geom_violin(aes(x = Group, y = WellExpressedPredFraction,
                  fill = Group, color = Group)) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) +
  ylab("Fraction expressed predictors") +
  coord_flip() +
  guides(color = F, fill = F, alpha = F) +
  theme_classic() + theme(text = element_text(size = 20),
                          legend.title = element_blank(),
                          axis.title.y = element_blank())


pdf(paste0("Plots/Human_Network/",net,"/Predictability/poorly_predicted/characterization.pdf"),
    width = 10)
ggarrange(plotlist = list(ggarrange(plotlist = list(ggarrange(plotlist = plots[1:2], ncol = 1, align = "hv"),
                          plots[[3]]), nrow = 1, common.legend = T, legend = 'none', labels = list("A)","B)")),
                          ggarrange(plotlist = plots[c(5,4)], nrow = 1, align = "h",
                                    common.legend = T, legend = "bottom", labels = list("C)","D)"))), nrow = 2,
          align = "hv", widths = c(5,4))
dev.off()


# GO term enrichment
go_bp_H <- GetGOEnrich(poorly_predicted_H, rownames(info), "BP", algorithm = "weight01", enrich_cutoff = 1)
go_bp_L <- GetGOEnrich(poorly_predicted_L, rownames(info), "BP", algorithm = "weight01", enrich_cutoff = 1)
go_mf_H <- GetGOEnrich(poorly_predicted_H, rownames(info), "MF", algorithm = "weight01", enrich_cutoff = 1)
go_mf_L <- GetGOEnrich(poorly_predicted_L, rownames(info), "MF", algorithm = "weight01", enrich_cutoff = 1)
go_cc_H <- GetGOEnrich(poorly_predicted_H, rownames(info), "CC", algorithm = "weight01", enrich_cutoff = 1)
go_cc_L <- GetGOEnrich(poorly_predicted_L, rownames(info), "CC", algorithm = "weight01", enrich_cutoff = 1)

go.plot.list <- list()
go.plot.list[[1]] <- PlotGOEnrich(go_bp_H, POORLY_PRED_H_COL,
                                  "Poorly Predicted\n(High expr.)")
go.plot.list[[2]] <- PlotGOEnrich(go_bp_L, POORLY_PRED_L_COL,
                                  "Poorly Predicted\n(Low expr.)")
go.plot.list[[3]] <- PlotGOEnrich(go_mf_H, POORLY_PRED_H_COL,
                                  "Poorly Predicted\n(High expr.)")
go.plot.list[[4]] <- PlotGOEnrich(go_mf_L, POORLY_PRED_L_COL,
                                  "Poorly Predicted\n(Low expr.)")
go.plot.list[[5]] <- PlotGOEnrich(go_cc_H, POORLY_PRED_H_COL,
                                  "Poorly Predicted\n(High expr.)")
go.plot.list[[6]] <- PlotGOEnrich(go_cc_L, POORLY_PRED_L_COL,
                                  "Poorly Predicted\n(Low expr.)")

pdf(paste0("Plots/Human_Network/",net,"/Predictability/poorly_predicted/GOmf.pdf"),
    height = 10, width = 14)
ggarrange(plotlist = go.plot.list[c(3,4)], ncol = 1, align = "v", heights = c(2,5))
dev.off()

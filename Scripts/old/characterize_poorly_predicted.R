POORLY_PRED_COL = "tomato2"
POORLY_PRED_H_COL = "tomato3"
POORLY_PRED_L_COL = "lightcoral"
COR_THRE = 0.2
LOW_MEAN_THRE = 2
HIGH_MEAN_THRE = 8

NET_FILE = "Data/Networks/Human/stabsel_network_Hs.rds"
net = "stabsel"

setwd("../")

# Load GTEx correlations
correlations <- readRDS(paste0("Outputs/Human_Network/",net,"/Predictability/Tissue/correlations_all.rds"))

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

# Load network
network_matrix <- as.matrix(readRDS(NET_FILE))
network_matrix <- network_matrix[rownames(GTEx_real),]

# # get mean and variance in training data
# train_data <- read.delim("/data/public/xwu2/CCTN_Scripts_and_DataSets/Data/CCTN_ConsideredCohortsData/CCLERSEQzcmb_GeneExpressionProfiles_24641x1443.txt")
# rownames(train_data) <- train_data[,1]
# train_data <- train_data[,-(1:3)]
# train_data <- train_data[names(GTEx_cor),]

# ----

plot.data <- data.frame("Genes" = rownames(GTEx_real), "Correlation" = correlations[,"CrossTissue"])
# Expression in test data
plot.data$GTExMeans <- rowMeans(GTEx_real)
plot(density(plot.data$GTExMeans))
abline(v = LOW_MEAN_THRE)
# Variance in test data
plot.data$GTExVars <- apply(GTEx_real, 1, var)
plot.data$Group <- c("Background","PoorlyPredicted")[as.numeric(plot.data$Correlation < COR_THRE)+1]

# # Expression in training data
# plot.data$TrainMeans <- rowMeans(train_data)
# plot(density(plot.data$TrainMeans))
# # Variance in training data
# plot.data$TrainRanges <- apply(train_data, 1, function(v) as.numeric(dist(range(v))))
# plot.data$TrainVars <- apply(train_data, 1, var)

library(reshape2)
library(ggplot2)
library(ggpubr)

pdf(paste0("Plots/Human_Network/",net,"/Predictability/poorly_predicted/mean_var_poorlypredicted.pdf"),
    height = 3, width = 5.5)

ggplot(plot.data) +
  geom_density(aes(x = GTExMeans), fill = "grey", color = "grey") +
  geom_density(data = subset(plot.data, Group == "PoorlyPredicted"),
               aes(x = GTExMeans), fill = POORLY_PRED_COL, color = POORLY_PRED_COL, alpha = .6) +
  xlab("Mean expression across GTEx samples") +
  # geom_vline(xintercept = 8) +
  theme_classic() + theme(text = element_text(size = 20),
    axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())

ggplot(plot.data) +
  geom_density(aes(x = GTExVars), fill = "grey", color = "grey") +
  geom_density(data = subset(plot.data, Group == "PoorlyPredicted"),
               aes(x = GTExVars), fill = POORLY_PRED_COL, color = POORLY_PRED_COL, alpha = .6) +
  xlab("Variance across GTEx samples") +
  theme_classic() + theme(text = element_text(size = 20),
    axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())

# ggplot(plot.data) +
#   geom_density(aes(x = TrainMeans), fill = "grey", color = "grey") +
#   geom_density(data = subset(plot.data, Group == "PoorlyPredicted"),
#                aes(x = TrainMeans), fill = POORLY_PRED_COL, color = POORLY_PRED_COL, alpha = .6) +
#   xlab("Mean expression across train samples") +
#   theme_classic() + theme(text = element_text(size = 20),
#     axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())
# 
# ggplot(plot.data) +
#   geom_density(aes(x = TrainVars), fill = "grey", color = "grey") +
#   geom_density(data = subset(plot.data, Group == "PoorlyPredicted"),
#                aes(x = TrainVars), fill = POORLY_PRED_COL, color = POORLY_PRED_COL, alpha = .6) +
#   xlab("Variance across train samples") +
#   theme_classic() + theme(text = element_text(size = 20),
#                           axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())

# ggplot(plot.data) +
#   geom_point(data = subset(plot.data, Group == "Background"),
#     aes(x = GTExMeans, y = GTExVars, color = Group), size = .5) +
#   geom_point(data = subset(plot.data, Group == "PoorlyPredicted"),
#     aes(x = GTExMeans, y = GTExVars, color = Group), size = .5) +
#   scale_color_manual(values = c("Background" = "grey", "PoorlyPredicted" = POORLY_PRED_COL)) +
#   xlab("Mean expression across GTEx samples") +
#   ylab("Expression var. across GTEx samples") +
#   theme_classic() + theme(text = element_text(size = 20))
dev.off()


# Number of predictors
plot.data$OriginalPredNumber <- rowSums(network_matrix[,] != 0)
network_matrix <- network_matrix[,intersect(colnames(network_matrix),rownames(GTEx_real))]
plot.data$QuantifiedPredNumber <- rowSums(network_matrix[,] != 0)

pdf(paste0("Plots/Human_Network/",net,"/Predictability/poorly_predicted/prednumb_poorlypredicted.pdf"),
    height = 4)
ggplot(subset(plot.data, !is.na(Correlation))) +
  geom_histogram(aes(x = OriginalPredNumber, fill = Group),
      bins = range(plot.data$OriginalPredNumber)[2]-range(plot.data$OriginalPredNumber)[1]) +
  scale_fill_manual(values = c("Background" = "grey", "PoorlyPredicted" = POORLY_PRED_COL)) +
  xlab("Number of predictors") + ylab("# genes") +
  theme_classic() + theme(text = element_text(size = 20), legend.title = element_blank())
ggplot(subset(plot.data,!is.na(Correlation))) +
  geom_histogram(aes(x = QuantifiedPredNumber, fill = Group),
      bins = range(plot.data$QuantifiedPredNumber)[2]-range(plot.data$QuantifiedPredNumber)[1]) +
  scale_fill_manual(values = c("Background" = "grey", "PoorlyPredicted" = POORLY_PRED_COL)) +
  xlab("Number of quantified predictors") + ylab("# genes") +
  theme_classic() + theme(text = element_text(size = 20), legend.title = element_blank())
dev.off()


# Predictor expression
plot.data$LowlyExpressedPredNumber <- sapply(rownames(plot.data), function(g){
  predictors <- names(which(network_matrix[g,] != 0))
  return(sum(plot.data$GTExMeans[plot.data$Genes %in% predictors] < 2))
})
plot.data$WellExpressedPredNumber <- sapply(rownames(plot.data), function(g){
  predictors <- names(which(network_matrix[g,] != 0))
  return(sum(plot.data$GTExMeans[plot.data$Genes %in% predictors] >= 2))
})
plot.data$WellExpressedPredFraction <- plot.data$WellExpressedPredNumber/plot.data$OriginalPredNumber
plot.data$LowlyExpressedPredFraction <- plot.data$LowlyExpressedPredNumber/plot.data$OriginalPredNumber

pdf(paste0("Plots/Human_Network/",net,"/Predictability/poorly_predicted/predexpr_poorlypredicted.pdf"),
    height = 4)

ggplot(subset(plot.data, !is.na(Correlation))) +
  geom_histogram(aes(x = WellExpressedPredNumber,
                     fill = Group),
                 bins = range(plot.data$WellExpressedPredNumber)[2]-
                   range(plot.data$WellExpressedPredNumber)[1]) +
  scale_fill_manual(values = c("Background" = "grey", "PoorlyPredicted" = POORLY_PRED_COL)) +
  xlab("Number of well expressed predictors") + ylab("# genes") +
  theme_classic() + theme(text = element_text(size = 20), legend.title = element_blank())

ggplot(subset(plot.data, !is.na(Correlation))) +
  geom_density(aes(x = WellExpressedPredFraction, fill = Group, color = Group, alpha = Group)) +
  scale_fill_manual(values = c("Background" = "grey", "PoorlyPredicted" = POORLY_PRED_COL)) +
  scale_color_manual(values = c("Background" = "grey", "PoorlyPredicted" = POORLY_PRED_COL)) +
  scale_alpha_manual(values = c("Background" = 1, "PoorlyPredicted" = .6)) +
  xlab("Fraction of well expressed predictors") + ylab("# genes") +
  guides(color = F, fill = F, alpha = F) +
  theme_classic() + theme(text = element_text(size = 20), legend.title = element_blank())

ggplot(plot.data) +
  geom_density(aes(x = LowlyExpressedPredFraction, fill = Group, color = Group, alpha = Group)) +
  scale_fill_manual(values = c("Background" = "grey", "PoorlyPredicted" = POORLY_PRED_COL)) +
  scale_color_manual(values = c("Background" = "grey", "PoorlyPredicted" = POORLY_PRED_COL)) +
  scale_alpha_manual(values = c("Background" = 1, "PoorlyPredicted" = .6)) +
  xlab("Fraction of lowly expressed predictors") + ylab("# genes") +
  guides(color = F, fill = F, alpha = F) +
  theme_classic() + theme(text = element_text(size = 20), legend.title = element_blank())

dev.off()

# most of the predictors are well expressed



# number of predictors in genes that are well captured vs not-well captured

plot.data$GroupDetail <- plot.data$Group
plot.data$GroupDetail[(plot.data$Group == "PoorlyPredicted") &
                        (plot.data$GTExMeans > HIGH_MEAN_THRE)] <- "PoorlyPredicted_H"
plot.data$GroupDetail[(plot.data$Group == "PoorlyPredicted") &
                        (plot.data$GTExMeans <= HIGH_MEAN_THRE)] <- "PoorlyPredicted_L"
plot.data$GroupDetail <- factor(plot.data$GroupDetail,
                                levels = c("Background","PoorlyPredicted_L","PoorlyPredicted_H"))
plots <- list()

plots[[1]] <- ggplot(subset(plot.data, !is.na(Correlation))) +
  geom_bar(aes(x = OriginalPredNumber, fill = GroupDetail)) +
  scale_fill_manual(values = c("Background" = "grey", "PoorlyPredicted_L" = POORLY_PRED_L_COL,
                               "PoorlyPredicted_H" = POORLY_PRED_H_COL)) +
  xlab("Number of predictors") + ylab("# genes") +
  theme_classic() + theme(text = element_text(size = 20), legend.title = element_blank())

# well captured but poorly predicted genes have a higher number of predictors... so issue must be smth else

plots[[2]] <- ggplot(subset(plot.data, !is.na(Correlation))) +
  geom_violin(aes(x = GroupDetail, y = WellExpressedPredFraction, fill = GroupDetail, color = Group)) +
  scale_fill_manual(values = c("Background" = "grey", "PoorlyPredicted_L" = POORLY_PRED_L_COL,
    "PoorlyPredicted_H" = POORLY_PRED_H_COL)) +
  scale_color_manual(values = c("Background" = "grey", "PoorlyPredicted_L" = POORLY_PRED_L_COL,
    "PoorlyPredicted_H" = POORLY_PRED_H_COL)) +
  scale_alpha_manual(values = c("Background" = 1, "PoorlyPredicted_L" = .6, "PoorlyPredicted_H" = 0.6)) +
  ylab("Fraction of well expressed predictors") +
  guides(color = F, fill = F, alpha = F) +
  theme_classic() + theme(text = element_text(size = 20), legend.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x = element_blank())

# well captured but poorly predicted genes also have most of their predictors expressed... so it has to be yet smth else

# Tissue-specificity of expression in test data
tissues <- readRDS("/data/public/adesous1/scDropImp/networkinference/analyses/paper/NetworkQuality_analysis/tissues.rds")
tissues <- droplevels(tissues)
tissue_means <- sapply(levels(tissues), function(t) rowMeans(GTEx_real[,as.character(tissues) == t]))
tissue_expressed <- rowSums(tissue_means > LOW_MEAN_THRE)
plot.data$TissuesExpressed <- tissue_expressed[rownames(plot.data)]
plots[[3]] <- ggplot(subset(plot.data, !is.na(Correlation))) +
  geom_bar(aes(x = TissuesExpressed, color = GroupDetail, fill = GroupDetail)) +
  scale_fill_manual(values = c("Background" = "grey", "PoorlyPredicted_L" = POORLY_PRED_L_COL,
    "PoorlyPredicted_H" = POORLY_PRED_H_COL)) +
  scale_color_manual(values = c("Background" = "grey", "PoorlyPredicted_L" = POORLY_PRED_L_COL,
    "PoorlyPredicted_H" = POORLY_PRED_H_COL)) +
  coord_cartesian(ylim = c(0,1000)) + xlab("# GTEx tissues where gene is expressed") + ylab("# genes") +
  guides(color = F, fill = F) + theme_classic() + theme(text = element_text(size = 20))

pdf(paste0("Plots/Human_Network/",net,"/Predictability/poorly_predicted/poorlypredicted_HvsL_preds_tissue.pdf"),
    width = 10)
ggarrange(plotlist = list(ggarrange(plotlist = plots[c(1,3)],
  ncol = 1, nrow = 2, common.legend = T, legend = "bottom"), plots[[2]]), ncol = 2, nrow = 1, widths = c(3,2))
dev.off()


# pdf("plots/poorlypredicted_HvsL_meanvartrain.pdf", height = 3, width = 5.5)
# ggplot(plot.data) +
#   geom_density(aes(x = TrainMeans, fill = GroupDetail, color = GroupDetail, alpha = GroupDetail)) +
#   xlab("Mean expression across train samples") + guides(fill = F, color = F, alpha = F) +
#   scale_fill_manual(values = c("Background" = "grey", "PoorlyPredicted_L" = POORLY_PRED_L_COL,
#     "PoorlyPredicted_H" = POORLY_PRED_H_COL)) +
#   scale_color_manual(values = c("Background" = "grey", "PoorlyPredicted_L" = POORLY_PRED_L_COL,
#     "PoorlyPredicted_H" = POORLY_PRED_H_COL)) +
#   scale_alpha_manual(values = c("Background" = 1, "PoorlyPredicted_L" = .5, "PoorlyPredicted_H" = .5)) +
#   theme_classic() + theme(text = element_text(size = 20),
#                           axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())
# 
# ggplot(plot.data) +
#   geom_density(aes(x = TrainVars, fill = GroupDetail, color = GroupDetail, alpha = GroupDetail)) +
#   scale_fill_manual(values = c("Background" = "grey", "PoorlyPredicted_L" = POORLY_PRED_L_COL,
#     "PoorlyPredicted_H" = POORLY_PRED_H_COL)) +
#   scale_color_manual(values = c("Background" = "grey", "PoorlyPredicted_L" = POORLY_PRED_L_COL,
#     "PoorlyPredicted_H" = POORLY_PRED_H_COL)) +
#   scale_alpha_manual(values = c("Background" = 1, "PoorlyPredicted_L" = .5, "PoorlyPredicted_H" = .5)) +
#   xlab("Variance across train samples") + guides(fill = F, color = F, alpha = F) +
#   theme_classic() + theme(text = element_text(size = 20),
#                           axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())
# dev.off()



pdf(paste0("Plots/Human_Network/",net,"/Predictability/poorly_predicted/poorlypredicted_HvsL_meanvarGTEx.pdf"),
    height = 3, width = 5.5)

ggplot(subset(plot.data, !is.na(Correlation))) +
  geom_density(aes(x = GTExMeans, fill = GroupDetail, color = GroupDetail, alpha = GroupDetail)) +
  xlab("Mean expression across GTEx samples") + guides(fill = F, color = F, alpha = F) +
  scale_fill_manual(values = c("Background" = "grey", "PoorlyPredicted_L" = POORLY_PRED_L_COL,
                               "PoorlyPredicted_H" = POORLY_PRED_H_COL)) +
  scale_color_manual(values = c("Background" = "grey", "PoorlyPredicted_L" = POORLY_PRED_L_COL,
                                "PoorlyPredicted_H" = POORLY_PRED_H_COL)) +
  scale_alpha_manual(values = c("Background" = 1, "PoorlyPredicted_L" = .5, "PoorlyPredicted_H" = .5)) +
  theme_classic() + theme(text = element_text(size = 20),
                          axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())

ggplot(subset(plot.data, !is.na(Correlation))) +
  geom_density(aes(x = GTExVars, fill = GroupDetail, color = GroupDetail, alpha = GroupDetail)) +
  scale_fill_manual(values = c("Background" = "grey", "PoorlyPredicted_L" = POORLY_PRED_L_COL,
                               "PoorlyPredicted_H" = POORLY_PRED_H_COL)) +
  scale_color_manual(values = c("Background" = "grey", "PoorlyPredicted_L" = POORLY_PRED_L_COL,
                                "PoorlyPredicted_H" = POORLY_PRED_H_COL)) +
  scale_alpha_manual(values = c("Background" = 1, "PoorlyPredicted_L" = .5, "PoorlyPredicted_H" = .5)) +
  xlab("Variance across GTEx samples") + guides(fill = F, color = F, alpha = F) +
  theme_classic() + theme(text = element_text(size = 20),
                          axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())
dev.off()

# GO enrichment highly expressed, ubiquitously expressed across tissues, poorly predicted

# background of ubiquitously expressed genes (>25/30)
background <- names(which(tissue_expressed[rownames(plot.data)] >=25))
# foreground of ubiquitously expressed and poorly predicted
foreground <- subset(plot.data, GroupDetail == "PoorlyPredicted_H")$Genes
source("Scripts/functions.R")
go_bp <- GetGOEnrich(foreground,background, "BP", algorithm = "weight01", enrich_cutoff = 1)
go_mf <- GetGOEnrich(foreground,background, "MF", algorithm = "weight01", enrich_cutoff = 1)
go_cc <- GetGOEnrich(foreground,background, "CC", algorithm = "weight01", enrich_cutoff = 1)

go.plot.list <- list()
go.plot.list[[1]] <- PlotGOEnrich(go_bp, POORLY_PRED_H_COL,
    "Poorly predicted, highly expressed vs ubiquitously expressed")
go.plot.list[[2]] <- PlotGOEnrich(go_mf, POORLY_PRED_H_COL,
    "Poorly predicted, highly expressed vs ubiquitously expressed")
go.plot.list[[3]] <- PlotGOEnrich(go_cc, POORLY_PRED_H_COL,
    "Poorly predicted, highly expressed vs ubiquitously expressed")

pdf(paste0("Plots/Human_Network/",net,"/Predictability/poorly_predicted/poorlypredicted_H_ubiquitous_GO.pdf"),
    width = 12, height = 10)
ggarrange(plotlist = go.plot.list[1:2], align = "hv", nrow = 2, heights = c(3.8,1.2))
dev.off()

pdf(paste0("Plots/Human_Network/",net,"/Predictability/poorly_predicted/poorlypredicted_H_ubiquitous_meanvarGTEx.pdf"),
    height = 3, width = 5.5)
ggplot() +
  geom_density(data = subset(plot.data, Genes %in% background),
               aes(x = GTExMeans), fill = "grey", color = "grey") +
  geom_density(data = subset(plot.data, Genes %in% foreground),
               aes(x = GTExMeans), fill = POORLY_PRED_H_COL, col = POORLY_PRED_L_COL, alpha = .5) +
  xlab("Mean expression across GTEx samples") +
  theme_classic() + theme(text = element_text(size = 20), axis.title.y = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggplot() +
  geom_density(data = subset(plot.data, Genes %in% background),
               aes(x = GTExVars), fill = "grey", color = "grey") +
  geom_density(data = subset(plot.data, Genes %in% foreground),
               aes(x = GTExVars), fill = POORLY_PRED_H_COL, col = POORLY_PRED_L_COL, alpha = .5) +
  xlab("Variance across GTEx samples") +
  theme_classic() + theme(text = element_text(size = 20), axis.title.y = element_blank(),
                          axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()


# background: genes highly expressed and ubiquitously expressed (What makes these genes hard to predict?)
background2 <- subset(plot.data, (GTExMeans > HIGH_MEAN_THRE) & (Genes %in% background))$Genes

go_bp_2 <- GetGOEnrich(foreground,background2, "BP", algorithm = "weight01", enrich_cutoff = 1)
go_mf_2 <- GetGOEnrich(foreground,background2, "MF", algorithm = "weight01", enrich_cutoff = 1)
go_cc_2 <- GetGOEnrich(foreground,background2, "CC", algorithm = "weight01", enrich_cutoff = 1)

PlotGOEnrich(go_bp_2, POORLY_PRED_H_COL, "Poorly predicted, highly expressed vs ubiquitously highly expressed")
PlotGOEnrich(go_mf_2, POORLY_PRED_H_COL, "Poorly predicted, highly expressed vs ubiquitously highly expressed")
PlotGOEnrich(go_cc_2, POORLY_PRED_H_COL, "Poorly predicted, highly expressed vs ubiquitously highly expressed")

# same crap


# background: all genes in network
# foreground: poorly predicted genes
background3 <- plot.data$Genes
foreground3 <- subset(plot.data, Group == "PoorlyPredicted")$Genes
go_bp_3 <- GetGOEnrich(foreground3,background3, "BP", algorithm = "weight01", enrich_cutoff = 1)
pdf(paste0("Plots/Human_Network/",net,"/Predictability/poorly_predicted/poorlypredicted_all_GObp.pdf"),
    width = 12, height = 8)
PlotGOEnrich(go_bp_3, POORLY_PRED_COL, "Poorly predicted vs all")
dev.off()

# foreground: poorly predicted lowly expressed genes
foreground_L <- subset(plot.data, GroupDetail == "PoorlyPredicted_L")$Genes
go_bp_L <- GetGOEnrich(foreground_L,background3, "BP", algorithm = "weight01", enrich_cutoff = 1)
pdf(paste0("Plots/Human_Network/",net,"/Predictability/poorly_predicted/poorlypredicted_L_all_GObp.pdf"),
    width = 11, height = 10)
PlotGOEnrich(go_bp_L, POORLY_PRED_L_COL, "Poorly predicted, lowly expressed vs all")
dev.off()

# foreground: poorly predicted highly expresssed genes
foreground_H <- subset(plot.data, GroupDetail == "PoorlyPredicted_H")$Genes
go_bp_H <- GetGOEnrich(foreground_H,background3, "BP", algorithm = "weight01", enrich_cutoff = 1)
pdf(paste0("Plots/Human_Network/",net,"/Predictability/poorly_predicted/poorlypredicted_H_all_GObp.pdf"),
    width = 11, height = 10)
PlotGOEnrich(go_bp_H, POORLY_PRED_H_COL, "Poorly predicted, highly expressed vs all")
dev.off()

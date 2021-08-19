# Plot predictability of genes not well predicted by stabsel network

NET_FILE = "Data/Networks/Human/stabsel_pcclasso_network_Hs_filtered.rds"
net = "stabsel_pcclasso"
POORLY_PRED_COL = "tomato2"

# Load all correlations

files <- list.files(paste0("Outputs/Human_Network/",net,"/Predictability/Tissue"),
                    pattern = "correlations.rds", full.names = T)
correlations <- sapply(files, readRDS)
colnames(correlations) <- sapply(files,
    function(x) strsplit(tail(strsplit(x, "/")[[1]], 1), "_")[[1]][1])
saveRDS(correlations, paste0("Outputs/Human_Network/",net,"/Predictability/Tissue/correlations_all.rds"))


# Load poorly predicted genes (according to net with partial correlations)

poorly_predicted <- readRDS("Outputs/Human_Network/stabsel/Predictability/Tissue/poorly_predicted_crosstissue.rds")


# Plot correlations of poorly and well predicted genes

plot.data <- melt(correlations)
plot.data$Group <- plot.data$Var1 %in% poorly_predicted
plot.data$Var2 <- gsub("(","\n(",plot.data$Var2,fixed=T)
plot.data$Var2 <- factor(as.character(plot.data$Var2),
                         levels = c("CrossTissue","Adipose-Subcutaneous","Artery-Tibial","Brain",
                                    "Esophagus-Mucosa","Lung","Muscle-Skeletal","Nerve-Tibial",
                                    "Skin-NotSunExposed\n(Suprapubic)","Skin-SunExposed\n(Lowerleg)",
                                    "Thyroid","WholeBlood"))

pdf(paste0("Plots/Human_Network/",net,
           "/Predictability/poorly_predicted_stabsel/CrossTissue_poorlypredicted_alltissues.pdf"),
    height = 8, width = 12)
ggplot(plot.data) +
  geom_violin(aes(y = value, x = Group, fill = Group), draw_quantiles = 0.5) +
  geom_hline(yintercept = 0.2, linetype = "dashed") +
  facet_wrap(~ Var2) + guides(fill = F) + xlab("") +
  ylab("Correlation coefficient") + ylim(c(-1,1)) +
  scale_fill_manual(values = c("TRUE" = POORLY_PRED_COL, "FALSE" = "darkgrey")) +
  scale_x_discrete(labels = c("TRUE" = "Poorly predicted", "FALSE" = "Others")) +
  theme(text = element_text(size = 20))
dev.off()

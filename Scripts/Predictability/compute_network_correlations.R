# Compute correlations between network predictions and original data

# PREDICTION_FOLDER = "Outputs/Human_Network/stabsel/Predictability/AgeTissue"
PREDICTION_FOLDER = "Outputs/Human_Network/stabsel/Predictability/Tissue"
# PREDICTION_FOLDER = "Outputs/Human_Network/stabsel_pcclasso/Predictability"
# PREDICTION_FOLDER = "Outputs/Human_Network/stabsel_pcclasso_filter01/Predictability/AgeTissue"
# PREDICTION_FOLDER = "Outputs/Human_Network/stabsel_randomized/Predictability/Tissue"

DATA_FOLDER = "GTEx_Networks/Tissue_Networks/Outputs"
# DATA_FOLDER = "GTEx_Networks/Age_Networks/Outputs"
# DATA_FOLDER = "GTEx_Networks/AgeTissue_Networks/Outputs"

files <- list.files(PREDICTION_FOLDER, pattern = "net_predictions.rds", full.names = T)
data_files <- list.files(DATA_FOLDER, pattern = "_sampled_centered_data.rds", full.names = T)

for(i in 1:length(files)){
  
  predicted <- readRDS(files[i])
  centered_data <- readRDS(data_files[i])
  
  # convert rownames to gene symbols
  conversion_table <- read.delim("Resources/ensembl_idversion_GTExDESeq2_symbolChrStart.txt")
  rownames(centered_data) <- sapply(rownames(centered_data),
                                    function(c) strsplit(c, split = "\\.")[[1]][1])
  rnames <- conversion_table[match(rownames(centered_data), conversion_table$ensembl_gene_id),"symbol"]
  centered_data <- centered_data[!is.na(rnames),]
  rownames(centered_data) <- rnames[!is.na(rnames)]
  rm(conversion_table,rnames); gc()
  
  centered_data <- centered_data[rownames(predicted),]
  
  if(sum(is.na(predicted)) > 0)
    break
  
  cor <- sapply(rownames(predicted), function(gene_name) cor.test(x = centered_data[gene_name,],
                                                                  y = predicted[gene_name,],
                                                                  method = "pearson")$estimate, simplify = T)
  names(cor) <- sapply(names(cor), function(x) strsplit(x,"\\.")[[1]][1])
  
  saveRDS(cor, gsub("net_predictions","correlations",files[i]))
}



# library(reshape2)
# library(ggplot2)
# library(viridis)
# 
# plot.data <- melt(cor)
# plot.data$Var2 <- factor(as.character(plot.data$Var2), levels = c("100","20","Original"))
# 
# # whole GTEx
# 
# pdf("netshuffling_GTExwhole_nointercept.pdf", width = 8, height = 5)
# ggplot() + 
#   geom_density(data = subset(plot.data, Var2 == "100"),
#                aes(x = value, color = Var2, fill = Var2)) +
#   geom_density(data = subset(plot.data, Var2 != "100"),
#                aes(x = value, color = Var2, fill = Var2), size = 3, alpha = 0.7) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
#   xlab("Correlation coefficient") +
#   scale_color_manual(values = c("Original" = "#1f73e0", "20" = "#ADD7F6", "100" = "darkgrey"),
#                      labels = c("Original" = "Network", "20" = "20% randomized network", "100" = "100% randomized network")) +
#   scale_fill_manual(values = c("Original" = "#1f73e0", "20" = "#ADD7F6", "100" = "darkgrey"),
#                     labels = c("Original" = "Network", "20" = "20% randomized network", "100" = "100% randomized network")) +
#   xlim(c(-1,1)) +
#   guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
#   theme_classic() + theme(text = element_text(size = 20),
#     legend.title = element_blank(), legend.position = "bottom",
#     axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
# dev.off()
# 
# 
# # excluding genes not expressed anywhere
# 
# GTEx_real <- GTEx_real[rownames(cor),]
# plot(density(rowMeans(GTEx_real)))
# expressed <- names(which(rowMeans(GTEx_real) > .5))
# cor_expressed <- cor[expressed,]
# saveRDS(cor_expressed, "cor_nointercept_whole_expressed.rds")
# 
# plot.data <- melt(cor_expressed)
# plot.data$Var2 <- factor(as.character(plot.data$Var2), levels = c("100","20","Original"))
# 
# pdf("netshuffling_GTExwhole_nointercept_expressed.pdf", width = 8, height = 5)
# ggplot() + 
#   geom_density(data = subset(plot.data, Var2 == "100"),
#                aes(x = value, color = Var2, fill = Var2)) +
#   geom_density(data = subset(plot.data, Var2 != "100"),
#                aes(x = value, color = Var2, fill = Var2), size = 3, alpha = 0.7) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
#   xlab("Correlation coefficient") +
#   scale_color_manual(values = c("Original" = "#1f73e0", "20" = "#ADD7F6", "100" = "darkgrey"),
#                      labels = c("Original" = "Network", "20" = "20% randomized network", "100" = "100% randomized network")) +
#   scale_fill_manual(values = c("Original" = "#1f73e0", "20" = "#ADD7F6", "100" = "darkgrey"),
#                     labels = c("Original" = "Network", "20" = "20% randomized network", "100" = "100% randomized network")) +
#   xlim(c(-1,1)) +
#   guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
#   theme_classic() + theme(text = element_text(size = 20),
#                           legend.title = element_blank(), legend.position = "bottom",
#                           axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
# dev.off()

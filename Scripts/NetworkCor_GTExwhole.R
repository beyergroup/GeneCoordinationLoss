setwd("/data/public/adesous1/scDropImp/networkinference/analyses/paper/NetworkQuality_analysis/GTEx_DESeq2/centered_allGTEx/correlation/")

files <- c("Original" = "/public-arc/xwu2/networks/infer_CCLERSEQ_01/predictions_GTEX_DESEQ2_centered_add1/Prediction_GTEx_ON_CCLERSEQ_01.expr.txt",
           "5" = "/public-arc/xwu2/networks/infer_CCLERSEQ_01/predictions_GTEX_DESEQ2_centered_add1_shuffle5/Prediction_GTEx_ON_CCLERSEQ_01.shuffle.per5.expr.txt",
           "20" = "/public-arc/xwu2/networks/infer_CCLERSEQ_01/predictions_GTEX_DESEQ2_centered_add1_shuffle20/Prediction_GTEx_ON_CCLERSEQ_01.shuffle.per20.expr.txt",
           "100" = "/public-arc/xwu2/networks/infer_CCLERSEQ_01/predictions_GTEX_DESEQ2_centered_add1_shuffle100/Prediction_GTEx_ON_CCLERSEQ_01.shuffle.per100.expr.txt")

library(data.table)
GTEx_real <- fread("/data/public/adesous1/scDropImp/networkinference/analyses/paper/NetworkQuality_analysis/GTEx_DESeq2/GTEx_DESeq2_GEinput.txt",
                   header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)
GTEx_real <- data.frame(GTEx_real,row.names = 1)
GTEx_real <- as.matrix(GTEx_real[,-(1:2)]) # dim 35488 17382
GTEx_real <- log2(GTEx_real + 1)


PredictionCor <- function(pred_file, real, add.intercept = F, add.centers = T){
  
  cat("Predicted expression from file:", pred_file, "\n")
  
  predicted <- fread(pred_file,
                     header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)
  predicted <- data.frame(predicted,row.names = 1)
  predicted <- as.matrix(predicted)
  
  if(add.intercept){
    predicted <- predicted[,-1]+predicted[,1]
  } else{
    predicted <- predicted[,-1]
  }
  
  # remove NAs
  predicted <- predicted[!(apply(predicted, 1, function(x) any(is.na(x)))),]
  
  # add centers
  if(add.centers){
    centers <- readRDS("../GTEx_DESeq2_centers.rds")
    names(centers) <- rownames(real)
    centers <- centers[rownames(predicted)]
    if(any(is.na(centers))){
      stop("NAs detected in gene expression centers\n")
    }
    predicted <- centers + predicted
  }

  # get common set of genes
  genes <- intersect(rownames(real), rownames(predicted))
  
  cor <- sapply(genes, function(gene_name) cor.test(x = real[gene_name,],
                                                    y = predicted[gene_name,],
                                                    method = "pearson")$estimate, simplify = T)
  gc()
  return(cor)
  
}

cor <- sapply(files[c(1,4)], PredictionCor, real = GTEx_real, add.intercept = F, simplify = T)
rownames(cor) <- sapply(rownames(cor), function(x) paste(head(strsplit(x, "\\.")[[1]],-1), collapse = "."))
saveRDS(cor, "cor_nointercept_whole.rds")

library(reshape2)
library(ggplot2)
library(viridis)

plot.data <- melt(cor)
plot.data$Var2 <- factor(as.character(plot.data$Var2), levels = c("100","20","Original"))

# whole GTEx

pdf("netshuffling_GTExwhole_nointercept.pdf", width = 8, height = 5)
ggplot() + 
  geom_density(data = subset(plot.data, Var2 == "100"),
               aes(x = value, color = Var2, fill = Var2)) +
  geom_density(data = subset(plot.data, Var2 != "100"),
               aes(x = value, color = Var2, fill = Var2), size = 3, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
  xlab("Correlation coefficient") +
  scale_color_manual(values = c("Original" = "#1f73e0", "20" = "#ADD7F6", "100" = "darkgrey"),
                     labels = c("Original" = "Network", "20" = "20% randomized network", "100" = "100% randomized network")) +
  scale_fill_manual(values = c("Original" = "#1f73e0", "20" = "#ADD7F6", "100" = "darkgrey"),
                    labels = c("Original" = "Network", "20" = "20% randomized network", "100" = "100% randomized network")) +
  xlim(c(-1,1)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_classic() + theme(text = element_text(size = 20),
    legend.title = element_blank(), legend.position = "bottom",
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()


# excluding genes not expressed anywhere

GTEx_real <- GTEx_real[rownames(cor),]
plot(density(rowMeans(GTEx_real)))
expressed <- names(which(rowMeans(GTEx_real) > .5))
cor_expressed <- cor[expressed,]
saveRDS(cor_expressed, "cor_nointercept_whole_expressed.rds")

plot.data <- melt(cor_expressed)
plot.data$Var2 <- factor(as.character(plot.data$Var2), levels = c("100","20","Original"))

# pdf("netshuffling_GTExwhole_nointercept_expressed.pdf", width = 8, height = 5)
pdf("Plots/Human_Network/stabsel/Predictability/netshuffling_GTExwhole_nointercept_expressed_adjusted.pdf", width = 6, height = 4)
ggplot() +
  geom_density(data = subset(plot.data, Var2 == "100"),
               aes(x = value, color = Var2, fill = Var2)) +
  # geom_density(data = subset(plot.data, Var2 != "100"),
  geom_density(data = subset(plot.data, Var2 == "Original"),
               aes(x = value, color = Var2, fill = Var2), size = 3, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
  xlab("Correlation coefficient") +
  # scale_color_manual(values = c("Original" = "#1f73e0", "20" = "#ADD7F6", "100" = "darkgrey"),
  #                    labels = c("Original" = "Network", "20" = "20% randomized network", "100" = "100% randomized network")) +
  # scale_fill_manual(values = c("Original" = "#1f73e0", "20" = "#ADD7F6", "100" = "darkgrey"),
  #                   labels = c("Original" = "Network", "20" = "20% randomized network", "100" = "100% randomized network")) +
  scale_color_manual(values = c("Original" = "#1f73e0", "100" = "darkgrey"),
                     labels = c("Original" = "True relationships", "100" = "Randomized relationships")) +
  scale_fill_manual(values = c("Original" = "#1f73e0", "100" = "darkgrey"),
                    labels = c("Original" = "True relationships", "100" = "Randomized relationships")) +
  xlim(c(-1,1)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_classic() + theme(text = element_text(size = 20),
                          legend.title = element_blank(), legend.position = "bottom",
                          axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()

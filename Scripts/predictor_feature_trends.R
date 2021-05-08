# Load GTEx correlations
GTEx_cor <- readRDS("/data/public/adesous1/scDropImp/networkinference/analyses/paper/NetworkQuality_analysis/GTEx_DESeq2/centered_allGTEx/correlation/cor_nointercept_whole_expressed.rds")
GTEx_cor <- GTEx_cor[,1]
GTEx_cor <- GTEx_cor[!is.na(GTEx_cor)]
names(GTEx_cor) <- sapply(names(GTEx_cor), function(x) strsplit(x, "\\.")[[1]][1])

# Load network
load("/data/public/adesous1/scDropImp/network/network.rds")
rm(O)
# Reduce to common targets
network_matrix <- network_matrix[intersect(rownames(network_matrix),names(GTEx_cor)),]

# Get initial number of predictors
original_pred_numbers <- rowSums(network_matrix != 0)

# Reduce to common predictors
network_matrix <- network_matrix[,intersect(colnames(network_matrix),names(GTEx_cor))]

# Get number of quantified predictors
quantified_pred_numbers <- rowSums(network_matrix != 0)

# Percentage of predictors that are quantified
quantified_pred_percent <- (quantified_pred_numbers/original_pred_numbers)*100

pdf("plots/prednumb_correlation_smoothscatter.pdf", height = 5, width = 5)
smoothScatter(x = original_pred_numbers, y = GTEx_cor[names(quantified_pred_numbers)],
              xlab = "Number of predictors", ylab = "Correlation in GTEx")
smoothScatter(x = quantified_pred_numbers, y = GTEx_cor[names(quantified_pred_numbers)],
              xlab = "Number of predictors captured in GTEx", ylab = "Correlation in GTEx")
smoothScatter(x = quantified_pred_percent, y = GTEx_cor[names(quantified_pred_percent)],
              xlab = "% predictors captured in GTEx", ylab = "Correlation in GTEx")
dev.off()


# -----

# Average expression of predictors

# get mean in GTEx
library(data.table)
GTEx_real <- fread("/data/public/adesous1/scDropImp/networkinference/analyses/paper/NetworkQuality_analysis/GTEx_DESeq2/GTEx_DESeq2_GEinput.txt",
                   header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)
GTEx_real <- data.frame(GTEx_real,row.names = 1)
GTEx_real <- as.matrix(GTEx_real[,-(1:2)]) # dim 35488 17382
GTEx_real <- log2(GTEx_real + 1)
GTEx_avgs <- rowMeans(GTEx_real[names(GTEx_cor),])
rm(GTEx_real); gc()

png("GTEx_avgs.png")
plot(density(GTEx_avgs)) # expr thre of 2
abline(v = 2)
dev.off()

THRE = 2

# Number of predictors well expressed
wellexpr_pred_numbers <- sapply(names(original_pred_numbers), function(g){
  predictors <- names(which(network_matrix[g,] != 0))
  return(sum(GTEx_avgs[predictors] > THRE))
})

# Percentage of predictors that are well expressed
wellexpr_pred_percent <- (wellexpr_pred_numbers/original_pred_numbers)*100

pdf("plots/predexpr_correlation_smoothscatter.pdf", height = 4, width = 4)
smoothScatter(x = wellexpr_pred_percent, y = GTEx_cor[names(wellexpr_pred_percent)],
              xlab = "% predictors that are well expressed in GTEx",
              ylab = "Correlation in GTEx")
dev.off()

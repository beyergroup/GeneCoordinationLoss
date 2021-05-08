# Load GTEx correlations
GTEx_cor <- readRDS("/data/public/adesous1/scDropImp/networkinference/analyses/paper/NetworkQuality_analysis/GTEx_DESeq2/centered_allGTEx/correlation/cor_nointercept_whole_expressed.rds")
GTEx_cor <- GTEx_cor[,1]
GTEx_cor <- GTEx_cor[!is.na(GTEx_cor)]

# Load TCGA correlations
TCGA_cor <- readRDS("/data/public/adesous1/scDropImp/networkinference/analyses/paper/NetworkQuality_analysis/TCGA/cor_pertissue/cor_nointercept.rds")
intersect <- lapply(TCGA_cor, rownames)
intersect <- do.call(c, intersect)
intersect <- table(intersect)
intersect <- names(intersect)[intersect == length(TCGA_cor)]
TCGA_cor <- lapply(TCGA_cor, function(m) m[intersect,1])
TCGA_cor <- do.call(cbind,TCGA_cor)
rm(intersect)
TCGA_cor <- TCGA_cor[rowSums(is.na(TCGA_cor)) != ncol(TCGA_cor),]


# reduce to common genes
rownames(TCGA_cor) <- sapply(rownames(TCGA_cor), function(x) strsplit(x, "\\.")[[1]][1])
names(GTEx_cor) <- sapply(names(GTEx_cor), function(x) strsplit(x, "\\.")[[1]][1])
intersect <- intersect(rownames(TCGA_cor),names(GTEx_cor))
TCGA_cor <- TCGA_cor[intersect,]
GTEx_cor <- GTEx_cor[intersect]


# get mean and variance in GTEx
library(data.table)
GTEx_real <- fread("/data/public/adesous1/scDropImp/networkinference/analyses/paper/NetworkQuality_analysis/GTEx_DESeq2/GTEx_DESeq2_GEinput.txt",
                   header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)
GTEx_real <- data.frame(GTEx_real,row.names = 1)
GTEx_real <- as.matrix(GTEx_real[,-(1:2)]) # dim 35488 17382
GTEx_real <- log2(GTEx_real + 1)
GTEx_avgs <- rowMeans(GTEx_real[intersect,])
GTEx_vars <- apply(GTEx_real[intersect,], 1, var)
rm(GTEx_real); gc()

# get mean and variance in training data
train_data <- read.delim("/data/public/xwu2/CCTN_Scripts_and_DataSets/Data/CCTN_ConsideredCohortsData/CCLERSEQzcmb_GeneExpressionProfiles_24641x1443.txt")
rownames(train_data) <- train_data[,1]
train_data <- train_data[,-(1:3)]
train_data <- train_data[intersect,]
train_avgs <- rowMeans(train_data)
train_ranges <- apply(train_data, 1, function(v) as.numeric(dist(range(v))))
train_vars <- apply(train_data, 1, var)
rm(train_data); gc()

# # save for GSEA
# cor_g <- sort(GTEx_cor, decreasing = TRUE)
# write.table(cbind(names(cor_g), cor_g), "GSEA/inputs/correlation.rnk", col.names = F, row.names = F,
#   quote = F, sep = "\t")

pdf("plots/correlation_trends.pdf", height = 4, width = 4)
smoothScatter(x = GTEx_avgs, y = GTEx_cor, xlab = "Mean expression GTEx", ylab = "Correlation in GTEx")
smoothScatter(x = log2(GTEx_vars), y = GTEx_cor, xlab = "Variance GTEx (log2)", ylab = "Correlation in GTEx")
smoothScatter(x = train_avgs, y = GTEx_cor, xlab = "Mean expression train", ylab = "Correlation in GTEx")
smoothScatter(x = log2(train_vars), y = GTEx_cor, xlab = "Variance train (log2)", ylab = "Correlation in GTEx")
dev.off()

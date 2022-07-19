# Compare top hits in Lung and Whole Blood using limma and correlation

library(limma, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(reshape2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggplot2, lib.loc = "Resources/Rlibs/R-4.0.3/")

tissue = "Lung"
tissue = "WholeBlood"

# Correlation hits
cor <- readRDS(paste0("Outputs/Human_Network/stabsel/Predictability/AgeTissue/",tissue,"_sampled_age_corP.rds"))
cor <- cor[order(cor[,"pval"], decreasing = F),]

# Limma hits
limma <- readRDS(paste0("Outputs/Human_Network/stabsel/Predictability/AgeTissue/",tissue,"_sampled_ageDP.rds"))
limma <- limma::topTable(fit = limma, coef = "AgeGroupOld", number = nrow(limma))

common <- intersect(rownames(limma), rownames(cor))

pdf(paste0("Plots/Human_Network/stabsel/Predictability/Age/tests/",tissue,"_logFC_cor.pdf"))
plot(x = limma[common,"logFC"], y = cor[common,"cor"],
     main = paste0(tissue), xlab = "logFC from distance linear model",
     ylab = "Correlation between age and pred-obs correlation")
abline(v = 0, col = "white")
dev.off()
# Results don't match :O


# Look at top hits from correlation (3 example genes)

genes <- rownames(cor)[1:3]
files <- list.files("Outputs/Human_Network/stabsel/Predictability/AgeTissue",
                    pattern = tissue, full.names = T)
files <- files[grep("sampled_correlations",files)]
correlations <- sapply(files, readRDS)
correlations <- correlations[genes,]
colnames(correlations) <- sapply(colnames(correlations),
                                 function(x) strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][2])

plot.data <- melt(correlations)
pdf(paste0("Plots/Human_Network/stabsel/Predictability/Age/tests/",tissue,"_top_cor_hits.pdf"))
ggplot(plot.data) +
  geom_point(aes(x = Var2, y = value, color = Var1)) +
  geom_line(aes(x = Var2, y = value, group = Var1, color = Var1)) +
  ylab("Correlation coefficient\n(predicted vs observed value)") +
  xlab("Age group") + ggtitle(paste0(tissue," top correlation hits")) +
  theme_classic() + theme(text = element_text(size = 20))
dev.off()

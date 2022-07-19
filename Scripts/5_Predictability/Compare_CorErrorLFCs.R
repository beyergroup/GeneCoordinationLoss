.libPaths("Resources/Rlibs/R-4.0.3/")
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
source("Scripts/functions.R")


error_files <- list.files("Outputs/5_Predictability",
                                   pattern = "_coefficients.rds",
                                   full.names = T)
cor_files <- list.files("Outputs/5_Predictability", pattern = "_corLFCs_YvsO.rds", full.names = T)
genecor_files <- list.files("Outputs/5_Predictability",
                            pattern = "_v7_sampled_centered_correlations.rds",
                            full.names = T)

well_predicted_genes <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")

tissues <- unique(sapply(error_files, function(x)
  strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1]))

plots <- list()

for(tissue in tissues){
  
  error <- readRDS(error_files[grepl(tissue,error_files,fixed =T)])
  cor <- readRDS(cor_files[grepl(tissue,cor_files,fixed =T)])
  
  error <- error[,"AgeSlope"]
  c <- readRDS(genecor_files[grepl(tissue,genecor_files,fixed =T)])
  c <- c[names(error)]
  error <- error*c
  
  common <- intersect(intersect(names(error),names(cor)),
                      well_predicted_genes[[tissue]])
  
  plot.data <- data.frame("Error slope" = error[common],
                          "Cor LFC" = cor[common])
  
  plots[[tissue]] <- ggplot(subset(plot.data, abs(Cor.LFC) < 5)) +
    geom_point(aes(y = Error.slope, x = Cor.LFC), alpha = .2) +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    theme_classic() + xlab("Cor.-based LFC") + ylab("Error-based LFC") +
    ggtitle(gsub("(","\n(",tissue,fixed=T)) + 
    theme(text = element_text(size = 20), title = element_text(size = 18))
  
}

pdf(paste0("Plots/5_Predictability/predictability_GTEx_subsets_tissue_LM_RMSE_corLFC.pdf"),
    height = 11, width = 16)
print(ggarrange(plotlist = plots, align = "hv"))
dev.off()

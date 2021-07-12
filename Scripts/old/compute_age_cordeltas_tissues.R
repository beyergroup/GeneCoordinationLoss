# Differential correlation with age in GTEx tissues

library(reshape2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggplot2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggpubr, lib.loc = "Resources/Rlibs/R-4.0.3/")

net = "stabsel"

setwd("../")

cor_files <- list.files(paste0("Outputs/Human_Network/",net,"/Predictability/AgeTissue"),
  pattern = "correlations.rds", full.names = T)
correlations <- sapply(cor_files, readRDS)

cor.data <- melt(correlations)
cor.data$Tissue <- sapply(as.character(cor.data$Var2),
  function(x) strsplit(tail(strsplit(x, "/")[[1]],1),"_")[[1]][1])
cor.data$Tissue <- sapply(cor.data$Tissue, function(x) gsub("(","\n(",x, fixed = T))
cor.data$Age <- sapply(as.character(cor.data$Var2),
  function(x) strsplit(tail(strsplit(x, "/")[[1]],1),"_")[[1]][2])
rm(correlations); gc()


# avg 50s and 60s - avg 20s and 30s
deltas <- list()
plots <- list()
for(tissue in unique(cor.data$Tissue)){
  
  y.data <- subset(cor.data, (Tissue == tissue) & (Age %in% c("20-29","30-39")))
  y.cors <- aggregate(y.data$value, by = list(y.data$Var1), FUN = function(x) mean(na.omit(x)))
  
  o.data <- subset(cor.data, (Tissue == tissue) & (Age %in% c("50-59","60-69")))
  o.cors <- aggregate(o.data$value, by = list(o.data$Var1), FUN = function(x) mean(na.omit(x)))
  
  rm(y.data, o.data); gc()
  cors <- merge.data.frame(y.cors, o.cors, by = "Group.1")
  cors$delta <- cors$x.y-cors$x.x
  
  deltas[[tissue]] <- cors$delta
  names(deltas[[tissue]]) <- cors$Group.1
  
  plots[[tissue]] <- eval(ggplot(cors) + geom_density(aes(x = delta)) +
    geom_vline(xintercept = 0, linetype = "dashed") + ggtitle(tissue) +
    xlim(c(-1,1)) + xlab("Age difference") + theme_classic() +
    theme(text = element_text(size = 20),
      legend.position = "bottom", legend.title = element_blank(),
      axis.title.y = element_blank(), title = element_text(size = 18)))
}

pdf(paste0("Plots/Human_Network/",net,"/Predictability/predictability_age_changes_tissues.pdf"),
    height = 9, width = 16)
ggarrange(plotlist = plots)
dev.off()

saveRDS(deltas, paste0("Outputs/Human_Network/",net,"/Predictability/AgeTissue/age_delta_tissues.rds"))


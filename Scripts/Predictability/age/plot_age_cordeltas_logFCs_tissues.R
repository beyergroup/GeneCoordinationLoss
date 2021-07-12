# Plot delta LFCs vs expression LFCs

library(limma, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(reshape2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggplot2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggpubr, lib.loc = "Resources/Rlibs/R-4.0.3/")

net = "stabsel"
# net = "stabsel_pcclasso_filter01"

predictability_files <- list.files("Outputs/Human_Network/stabsel/Predictability/AgeTissue",
  pattern = "_sampled_ageDP.rds", full.names = T)
expression_files <- list.files("Outputs/GTEx/AgeDE",
  pattern = "_ageDE.rds", full.names = T)

tissues <- unique(sapply(predictability_files, function(x)
  strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1]))

plots <- list()
for(tissue in tissues){
  
  ageDE <- readRDS(expression_files[grepl(tissue,expression_files,fixed =T)])
  ageDE <- limma::topTable(fit = ageDE, coef = "AgeGroupOld", number = nrow(ageDE))
  
  # convert rownames to gene symbols
  conversion_table <- read.delim("Resources/ensembl_idversion_GTExDESeq2_symbolChrStart.txt")
  rownames(ageDE) <- sapply(rownames(ageDE),
    function(c) strsplit(c, split = "\\.")[[1]][1])
  rnames <- conversion_table[match(rownames(ageDE), conversion_table$ensembl_gene_id),"symbol"]
  duplicated <- na.omit(rnames)[duplicated(na.omit(rnames))]
  ageDE <- ageDE[(!is.na(rnames)) & (!(rnames %in% duplicated)),]
  rownames(ageDE) <- rnames[(!is.na(rnames)) & (!(rnames %in% duplicated))]
  rm(conversion_table,rnames); gc()
  
  ageDP <- readRDS(predictability_files[grepl(tissue,predictability_files,fixed = T)])
  ageDP <- limma::topTable(fit = ageDP, coef = "AgeGroupOld", number = nrow(ageDP))
  
  common <- intersect(rownames(ageDP),rownames(ageDE))
  
  plot.data <- data.frame("DE_LFC" = ageDE[common,"logFC"], "DE_pval" = ageDE[common,"adj.P.Val"],
                          "DP_LFC" = ageDP[common,"logFC"], "DP_pval" = ageDP[common,"adj.P.Val"])
  
  plots[[tissue]] <- ggplot() +
    geom_point(data = subset(plot.data, DE_pval > 0.05 & DP_pval > 0.05),
               aes(x = DE_LFC, y = DP_LFC), size = .5, color = "grey") +
    geom_point(data = subset(plot.data, DE_pval < 0.05),
               aes(x = DE_LFC, y = DP_LFC), color = "darkgrey") +
    geom_point(data = subset(plot.data, DP_pval < 0.05),
               aes(x = DE_LFC, y = DP_LFC), color = "orangered2") +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    ylim(c(-max(abs(plot.data$DP_LFC)),max(abs(plot.data$DP_LFC)))) +
    xlim(c(-max(abs(plot.data$DE_LFC)),max(abs(plot.data$DE_LFC)))) +
    xlab("Expression fold change") + ylab("Deviation fold change") +
    guides(alpha = F, color = F) + theme_classic() +
    ggtitle(gsub("(","\n(",tissue,fixed=T)) + 
    theme(text = element_text(size = 20), title = element_text(size = 18))
}

pdf(paste0("Plots/Human_Network/",net,"/Predictability/predictability_GTEx_subsets_tissue_age_DE.pdf"),
  height = 10, width = 16)
ggarrange(plotlist = plots)
dev.off()

# Boxplots of LFCs for individual TADs


library(ggplot2)
source("Scripts/functions.R")


tissues <- list("Artery-Tibial" = "Aorta", "Lung" = c("A549","Lung"))

# Well-predicted genes
well_predicted_genes <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")


for(tissue in names(tissues)){
  
  files <- list.files("Outputs/6_GenomicLocation",
                      pattern = paste(tissues[[tissue]],collapse="|"),
                      full.names = TRUE)
  files <- files[grep("TADs_genes.rds",files)]
  
  for(file in files){
    
    # Load mapping of TADs to genes
    TAD2gene <- ReadRDS(file)
    
    # Load tissue age LFCs
    ageLFCs <- ReadRDS(paste0("Outputs/5_Predictability/",tissue,"_trans_corLFCs_YvsO_adj_spearman.rds"))
    
    # Make boxplot
    plot.data <- subset(TAD2gene,
                        hgnc_symbol %in% intersect(names(ageLFCs),
                                                   well_predicted_genes[[tissue]]),
                        select = c("TAD_ID","hgnc_symbol","Chromosome"))
    plot.data$LFCs <- ageLFCs[plot.data$hgnc_symbol]
    plot.data$Chromosome <- factor(as.character(plot.data$Chromosome),
                                   levels = c(as.character(1:22),"X"))
    plot.data$TAD_ID <- factor(as.character(plot.data$TAD_ID),
                               levels = paste0("TAD",
                                               1:max(as.numeric(gsub("TAD","",
                                                                     TAD2gene$TAD_ID)))))
    plot.data <- subset(plot.data,
                        TAD_ID %in% names(which(table(plot.data$TAD_ID) > 2)))
    
    # Test for difference from 0
    pvals <- c()
    for(i in unique(plot.data$TAD_ID)){
      pvals <- c(pvals, wilcox.test(subset(plot.data, TAD_ID == i)$LFCs,
                                    alternative = "two.sided")$p.val)
    }
    pvals <- p.adjust(pvals, method = "fdr")
    names(pvals) <- unique(plot.data$TAD_ID)
    
    if(any(pvals < 0.1)){
      message("Significantly altered TADs:", names(which(pvals < 0.1)))
    }
    
    pdf(paste0("Plots/6_GenomicLocation/",
               gsub("_genes.rds","",tail(strsplit(file,"/")[[1]],1)),
               "_TAD_trans_LFC_boxplots.pdf"), height = 28, width = 12)
    print(ggplot(plot.data) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            geom_boxplot(aes(x = TAD_ID, y = LFCs),
                         size = .1, fill = "transparent", outlier.size = 0) +
            geom_jitter(aes(x = TAD_ID, y = LFCs),
                        size = .5, width = .1) +
            facet_wrap(~ Chromosome, ncol = 1, scales = "free_x",
                       labeller = label_both) +
            xlab("TADs") + ylab("Age error LFC") +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank()))
    dev.off()
  }
}

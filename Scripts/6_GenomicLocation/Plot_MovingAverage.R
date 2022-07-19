# Plot moving average of LFCs for all tissues

library(ggplot2)
source("Scripts/functions.R")


WINDOW = "1"


moving_averages <- ReadRDS("Outputs/6_GenomicLocation/moving_averages_trans.rds")

moving_averages <- do.call(rbind, moving_averages)


pdf(paste0("Plots/6_GenomicLocation/MA",WINDOW,"_trans_alltissues.pdf"),
    width = 20, height = 20)

for(chr in c(as.character(1:22),"X")){
  
  plot.data <- subset(moving_averages, Chromosome == chr)
  
  print(ggplot(plot.data) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
          geom_line(aes_string(x = "Position", y = paste0("MA",WINDOW))) +
      facet_wrap(~ Tissue, ncol = 1) +
      xlab("Genomic position (bp)") +
      ylab("Age-related error LFC") +
      ggtitle(paste("Chromosome",chr)) +
      scale_x_continuous(expand = c(0,0)) +
      theme_bw() + theme(text = element_text(size = 20)))
}

dev.off()

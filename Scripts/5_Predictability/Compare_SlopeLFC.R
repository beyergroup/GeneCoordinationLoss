source("Scripts/functions.R")
library(ggplot2)

slope_files <- list.files("Outputs/5_Predictability/tmp",
                          pattern = "ageslope_0.75", full.names = T)
slopes <- sapply(slope_files, ReadRDS, simplify = F)
lapply(slopes, function(m) sum(m[,"FDR"] < 0.1)) # 1 gene max in some tissues
names(slopes) <- sapply(names(slopes),
                        function(x) strsplit(tail(strsplit(x, "/")[[1]],1),
                                             "_")[[1]][1])
rm(slope_files); gc()

# compare slopes with LFCs
lfc_files <- list.files("Outputs/5_Predictability",
                        pattern = "YvsO", full.names = T)
lfcs <- sapply(lfc_files, ReadRDS, simplify = F)
names(lfcs) <- sapply(names(lfcs),
                      function(x) strsplit(tail(strsplit(x, "/")[[1]],1),
                                           "_")[[1]][1])
rm(lfc_files); gc()


plot.data <- data.frame()
for(tissue in names(slopes)){
  
  genes <- intersect(names(lfcs[[tissue]]), rownames(slopes[[tissue]]))
  plot.data <- rbind.data.frame(plot.data,
                                data.frame("Slope" = slopes[[tissue]][genes,
                                                                      "AdjustedSlope"],
                                           "LFC" = lfcs[[tissue]][genes],
                                           "Tissue" = tissue))
}

pdf("Plots/5_Predictability/tmp/slopes_LFCs.pdf", width = 11, height = 15)
ggplot(plot.data) +
  geom_point(aes(x = LFC, y = Slope), alpha = .2) +
  facet_wrap(~ Tissue, ncol = 3) +
  theme(text = element_text(size = 20))
dev.off()

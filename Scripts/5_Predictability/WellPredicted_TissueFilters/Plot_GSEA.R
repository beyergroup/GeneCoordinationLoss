.libPaths("Resources/Rlibs/R-4.3.1")
library(ggplot2)
source("Scripts/functions.R")

files <- list.files("Outputs/5_Predictability/WellPredicted_TissueFilters/GSEA/Out",
                    pattern = "gsea_report", recursive = TRUE,
                    full.names = TRUE)
files <- files[grep(".tsv",files)]

type = "GOmf"

files <- files[grep(type,files)]

tissues <- unique(sapply(files,
                         function(x) strsplit(strsplit(x, split = "/")[[1]][7],
                                              "\\.")[[1]][1]))

pdf(paste0("Plots/5_Predictability/WellPredicted_TissueFilters/GSEA_",type,".pdf"),
    width = 10, heigh = 10)

for(tissue in tissues){
  
  enrichments <- sapply(files[grep(tissue,files)],
                        function(x) read.delim(x, sep = "\t")[,c("NAME","NES","FDR.q.val")],
                        simplify = FALSE)
  enrichments <- do.call(rbind.data.frame, enrichments)
  enrichments <- subset(enrichments, FDR.q.val < 0.05)
  enrichments <- enrichments[order(enrichments$NES, decreasing = T),]
  
  if(nrow(enrichments) > 20){
    enrichments <- rbind.data.frame(head(enrichments,10),
                                    tail(enrichments,10))
  }
  enrichments$NAME <- factor(as.character(enrichments$NAME),
                             levels = enrichments$NAME)
  enrichments$NES <- as.numeric(enrichments$NES)
  
  print(ggplot(enrichments) +
    geom_bar(aes(x = NAME, y = NES, alpha = -log10(as.numeric(FDR.q.val)),
                 fill = as.numeric(NES)),
             stat = "identity") + guides(alpha = "none") +
    geom_vline(xintercept = 0) +
    scale_fill_gradient2(low ="tomato2", mid = "white", high = "#354f52", 
                         midpoint = 0, name = "NES") +
    ggtitle(tissue) + coord_flip() + theme_classic() +
    theme(text = element_text(size = 20)))
  
}

dev.off()

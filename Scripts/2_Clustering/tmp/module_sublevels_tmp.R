source("Scripts/functions.R")
source("Scripts/2_Clustering/params.R")
library(ggpubr)

membership <- readRDS("Outputs/2_Clustering/stabsel_filtered_largestCC_Adjacency_undirected_weightssum_rownormalized_membership_reclustering_redecomposing_250_adjusted.rds")
head(sort(table(membership), decreasing = T))

membership_lvl1 <- sapply(membership, function(x) strsplit(x, "\\.")[[1]][1])
membership_lvl2 <- sapply(membership, function(x) paste(strsplit(x, "\\.")[[1]][1:2], collapse = "."))
membership_lvl3 <- sapply(membership, function(x) paste(strsplit(x, "\\.")[[1]][1:3], collapse = "."))
membership_lvl4 <- sapply(membership, function(x) paste(strsplit(x, "\\.")[[1]][1:4], collapse = "."))
membership_lvl5 <- sapply(membership, function(x) paste(strsplit(x, "\\.")[[1]][1:5], collapse = "."))
membership_lvl6 <- sapply(membership, function(x) paste(strsplit(x, "\\.")[[1]][1:6], collapse = "."))
membership_lvl7 <- sapply(membership, function(x) paste(strsplit(x, "\\.")[[1]][1:7], collapse = "."))
membership_lvl8 <- sapply(membership, function(x) paste(strsplit(x, "\\.")[[1]][1:8], collapse = "."))
membership_lvl9 <- sapply(membership, function(x) paste(strsplit(x, "\\.")[[1]][1:9], collapse = "."))
membership_lvl10 <- sapply(membership, function(x) paste(strsplit(x, "\\.")[[1]][1:10], collapse = "."))
membership_lvl11 <- sapply(membership, function(x) paste(strsplit(x, "\\.")[[1]][1:11], collapse = "."))

load("/data/public/adesous1/GeneCorrelation/Outputs/2_Clustering/stabsel_filtered_largestCC_Adjacency_undirected_weightssum_rownormalized_membership_reclustering_redecomposing_250_adjusted_GO_weight01.RData")

modules = names(which(sort(table(membership),decreasing = T) >= 10))
WriteRDS(modules, "Outputs/2_Clustering/final_modules.rds")
modules <- c(modules, "1.2", "1.2.3","1.2.3.1","1.2.3.1.1.2","1.2.2","1.2.4","1.2.5",
             "1.1","1.1.2","1.1.2.1.1",
             "5.1.1", "8.1", "13.1")

for(MODULE in modules){
  
  if(!(MODULE %in% names(GObp.list))){
    
    tmp_membership <- 
    genes <- names(which(membership_lvl6 == MODULE))
    
    GObp.list[[MODULE]] <- GetGOEnrich(genes, names(membership), "BP", enrich_cutoff = 1, algorithm = "weight01")
    GOmf.list[[MODULE]] <- GetGOEnrich(genes, names(membership), "MF", enrich_cutoff = 1, algorithm = "weight01")
    GOcc.list[[MODULE]] <- GetGOEnrich(genes, names(membership), "CC", enrich_cutoff = 1, algorithm = "weight01")
  }
  # next module: 114
  MODULE = modules[i]
  if(sum(nrow(GObp.list[[MODULE]]),nrow(GOmf.list[[MODULE]]),nrow(GOcc.list[[MODULE]])) == 0)
    next
  
  plot.list <- list()
  plot.list[[1]] <- PlotGOEnrich(GObp.list[[MODULE]], col = COLORS[MODULE],
                                 title = paste0("GObp in module ",MODULE))
  plot.list[[2]] <- PlotGOEnrich(GOmf.list[[MODULE]], col = COLORS[MODULE],
                                 title = paste0("GOcc in module ",MODULE))
  plot.list[[3]] <- PlotGOEnrich(GOcc.list[[MODULE]], col = COLORS[MODULE],
                                 title = paste0("GOmf in module ",MODULE))
  
  heights <- c(nrow(GObp.list[[MODULE]]), nrow(GOmf.list[[MODULE]]), nrow(GOcc.list[[MODULE]])) + 3.5
  
  pdf(paste0("Plots/2_Clustering/tmp/GO_weight01_module_",MODULE,".pdf"),
      width = PDF_DIMS[MODULE,"Width"], height = PDF_DIMS[MODULE,"Height"])
  ggarrange(plotlist = plot.list, ncol = 1, nrow = length(plot.list),
            heights = heights, align = "hv")
  dev.off()
}

sizes <- as.numeric(table(membership))
names(sizes) <- names(table(membership))
sizes <- sizes[sizes >= 10]
sizes <- sort(sizes, decreasing = T)

pdf(paste0("Plots/2_Clustering/tmp/iterative_reclustering_redecomposing_module_sizes.pdf"),
    height = 5, width = 8)
print(barplot(sizes, xlab = "Module ID", ylab = "Module size", col = "snow3",
              main = "Iterative reclustering", sub = "Max size = 250"))
dev.off()



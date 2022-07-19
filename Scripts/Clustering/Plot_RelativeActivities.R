# Plot relative activities

.libPaths("Resources/Rlibs/R-4.0.3/")
library(pheatmap)
library(RColorBrewer)

DEND_FILE = "Outputs/Human_Network/stabsel/Modules/old/Adjacency_weightnone_pruned_dendrogram.rds"
INDIR = "Outputs/Human_Network/stabsel/Modules/Module_Activity/"
# INDIR = "Outputs/Human_Network/stabsel/Modules/Module_Activity/LocalWeightedActivity/"
# INDIR = "Outputs/Human_Network/stabsel/Modules/Module_Activity/GlobalSmoothedActivity/"
# OUTDIR = "Plots/Human_Network/stabsel/Modules/Module_Activity/GlobalSmoothedActivity/"
# OUTDIR = "Plots/Human_Network/stabsel/Modules/Module_Activity/LocalWeightedActivity/"
OUTDIR = "Plots/Human_Network/stabsel/Modules/Module_Activity/"

dir.create(OUTDIR)

# Read in dendrogram of clusters ----------------------------------------------

dend <- readRDS(DEND_FILE)

# Read in files of averaged smoothed values -----------------------------------

files <- list.files(INDIR, pattern = ".rds")

for(file in files){
  
  activity <- readRDS(paste0(INDIR,file))
  activity <- activity[,-grep("CrossTissue",colnames(activity))]
  plot.data <- as.matrix(activity)
  
  pheatmap(t(plot.data), scale = "none", border_color = NA,
           color = colorRampPalette(c("blue","white","red"))(n = 501),
           breaks = seq(from = -max(abs(plot.data)), to = max(abs(plot.data)), length.out = 501),
           cellwidth = 10, cellheight = 10, show_colnames = T,
           cluster_cols = as.hclust(dend),
           main = "Relative module activity in GTEx tissues",
           filename = paste0(OUTDIR, gsub(".rds",".pdf",file)))
}

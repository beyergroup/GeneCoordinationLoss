# Evaluate network modules

WDIR = "/data/public/adesous1/GeneCorrelation/"
setwd(WDIR)

.libPaths("Resources/Rlibs/R-4.0.3/")
library(igraph)
library(pheatmap)
library(RColorBrewer)
library(topGO)
source("Scripts/functions.R")

net = "stabsel"
type = "network_largest_cc"
weights = "none"
matrix = "Laplacian"


# read in decisions
decisions <- readRDS(paste0("Outputs/Human_Network/",net,"/Modules/",matrix,
                            "_weight",weights,"_height_eigen_decision.rds"))

# read in corresponding hc object
hc <- readRDS(paste0("Outputs/Human_Network/",net,"/Modules/",matrix,"_weight",
                     weights,"_complete_link_clustering_",
                     as.character(decisions$Eigenvectors),"_evectors.rds"))

# cut at chosen height
membership <- cutree(hc, h = decisions$Height)

# keep dendrogram for plotting
dend <- as.dendrogram(hc)
dend_cut <- cut(dend, h = decisions$Height)
plot(dend_cut$upper)

# cluster size distribution
clusters <- table(membership)
barplot(sort(clusters, decreasing = TRUE))
min(clusters)
max(clusters)
pdf(paste0("Plots/Human_Network/",net,"/Modules/",matrix,"_weight",weights,
           "_complete_link_cluster_size.pdf"))
barplot(sort(clusters[clusters > 2], decreasing = TRUE),
        ylab = "Cluster size ", xlab = "Clusters", log = "y", names.arg = "",
        main = paste0(sum(clusters > 2)," out of ", length(clusters)," clusters"),
        sub = "Cluster sizes < 3 excluded")
plot(ecdf(sort(clusters)))
dev.off()


modules <- names(table(membership)[table(membership) >= 10])

# cluster expression and variance in tissues

tissue_files <- list.files("GTEx_Networks/Tissue_Networks/Outputs",
                           pattern = "sampled_data.rds", full.names = T)

gtex_norm <- readRDS("Data/GTEx/DESeq2_normalized_gtex.rds")


abs_expression <- list()
variance <- list()
within_correlation <- list()

for(file in tissue_files){
  
  # read in tissue-specific GTEx subset (for sample identification)
  centered_gtex <- readRDS(file)
  centered_gtex <- DataENSGToSymbol(centered_gtex, remove_dup = T)
  samples <- colnames(centered_gtex)
  
  # subset big GTEx
  gtex <- gtex_norm[,samples]
  rm(samples); gc()
  # convert to gene symbol
  gtex <- DataENSGToSymbol(gtex, remove_dup = T)
  gtex <- log2(1+gtex)
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  abs_expression[[tissue]] <- c()
  variance[[tissue]] <- c()
  within_correlation[[tissue]] <- c()
  
  for(module in modules){
    genes <- names(which(membership == module))
    abs_expression[[tissue]] <- c(abs_expression[[tissue]],
                                  mean(gtex[intersect(genes,rownames(gtex)),]))
    variance[[tissue]] <- c(variance[[tissue]],
                            mean(apply(gtex[intersect(genes,rownames(gtex)),],1,var)))
    
    cor_mat <- cor(t(centered_gtex[intersect(genes,rownames(centered_gtex)),]),
                   method = "pearson")
    
    within_correlation[[tissue]] <- c(within_correlation[[tissue]],
                                      mean(abs(cor_mat[upper.tri(cor_mat)])))
  }
  names(abs_expression[[tissue]]) <- modules
  names(variance[[tissue]]) <- modules
  names(within_correlation[[tissue]]) <- modules
}

abs_expression <- do.call(cbind, abs_expression)
variance <- do.call(cbind, variance)
within_correlation <- do.call(cbind, within_correlation)

ann <- data.frame("Module size" = as.numeric(table(membership)[modules]))
rownames(ann) <- modules

# expression heatmap
pheatmap(t(abs_expression), scale = "none", border_color = NA,
         color = colorRampPalette(c("white","red"))(n = 501),
         cellwidth = 10, cellheight = 10, cluster_cols = F,
         annotation_col = ann, annotation_colors = list("Module.size" = c("white","black")),
         main = paste0("Absolute expression (average per module)\n",
                       decisions$Eigenvectors," ",matrix," eigenvectors, ",c("none" = "unweighted",
                                                                             "abs" = "|weights|")[weights]),
         filename = paste0("Plots/Human_Network/",net,"/Modules/Expression/absexpression_",
                           matrix,"_weight",weights,".pdf"))

# variance heatmap
pheatmap(t(variance), scale = "none", border_color = NA,
         color = colorRampPalette(c("white","red"))(n = 501),
         cellwidth = 10, cellheight = 10, cluster_cols = F,
         annotation_col = ann, annotation_colors = list("Module.size" = c("white","black")),
         main = paste0("Expression variance (average per module)\n",
                       decisions$Eigenvectors," ",matrix," eigenvectors, ",c("none" = "unweighted",
                                                                             "abs" = "|weights|")[weights]),
         filename = paste0("Plots/Human_Network/",net,"/Modules/Expression/var_",
                           matrix,"_weight",weights,".pdf"))

# within-correlation heatmap
pheatmap(t(within_correlation), scale = "none", border_color = NA,
         # color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdWhBu")))(n = 501),
         colorRampPalette(c("blue","white","red"))(n = 501),
         breaks = seq(from = -1, to = 1, length.out = 501),
         cellwidth = 10, cellheight = 10, cluster_cols = F,
         annotation_col = ann, annotation_colors = list("Module.size" = c("white","black")),
         main = paste0("Within-module correlation (average per module)\n",
                       decisions$Eigenvectors," ",matrix," eigenvectors, ",c("none" = "unweighted",
                                                                             "abs" = "|weights|")[weights]),
         filename = paste0("Plots/Human_Network/",net,"/Modules/Expression/withincorrelation_",
                           matrix,"_weight",weights,".pdf"))



# predictability

predictability_files <- list.files(paste0("Outputs/Human_Network/",net,"/Predictability/Tissue"),
                                   pattern = "sampled_correlations.rds", full.names = T)

predictability <- list()

for(file in predictability_files){
  
  pred <- readRDS(file)
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  predictability[[tissue]] <- c()
  
  for(module in modules){
    genes <- names(which(membership == module))
    predictability[[tissue]] <- c(predictability[[tissue]],
                                  mean(pred[intersect(genes,rownames(gtex))], na.rm = T))
  }
  names(predictability[[tissue]]) <- modules
}

predictability <- do.call(cbind, predictability)

ann <- data.frame("Module size" = as.numeric(table(membership)[modules]))
rownames(ann) <- modules

# predictability heatmap
pheatmap(t(predictability), scale = "none", border_color = NA,
         # color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(n = 501),
         colorRampPalette(c("blue","white","red"))(n = 501),
         breaks = seq(from = -1, to = 1, length.out = 501),
         cellwidth = 10, cellheight = 10, cluster_cols = F,
         annotation_col = ann, annotation_colors = list("Module.size" = c("white","black")),
         main = paste0("Predictability (average per module)\n",
                       decisions$Eigenvectors," ",matrix," eigenvectors, ",c("none" = "unweighted",
                                                                             "abs" = "|weights|")[weights]),
         filename = paste0("Plots/Human_Network/",net,"/Modules/Expression/predictability_",
                           matrix,"_weight",weights,".pdf"))


# GO term enrichment

GObp.list <- list()
GOmf.list <- list()
GOcc.list <- list()

for(module in modules){
  
  genes <- names(which(membership == module))
  
  GObp.list[[module]] <- GetGOEnrich(genes, names(membership), "BP", enrich_cutoff = 1, algorithm = "elim")
  GOmf.list[[module]] <- GetGOEnrich(genes, names(membership), "MF", enrich_cutoff = 1, algorithm = "elim")
  GOcc.list[[module]] <- GetGOEnrich(genes, names(membership), "CC", enrich_cutoff = 1, algorithm = "elim")
  
}
save(GObp.list,GOmf.list,GOcc.list, file = paste0("Outputs/Human_Network/",net,"/Modules/GO_enrichment/",matrix,
                                                  "_weight",weights,"_elim.RData"))

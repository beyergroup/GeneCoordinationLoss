# Read in GSEA results

args <- commandArgs(trailingOnly = TRUE)
GSEA_DIR = args[1]
TYPE = args[2]
PLOT_DIR = args[3]

library(readr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(gtools)
source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")

type = "Hallmarks"


files <- list.files(paste0("Outputs/5_Predictability/Age/Age_MaxSubset/SmoothedSlopes/GSEA/Outputs/",type),
                    recursive = T, pattern = "gsea_report_for_na_")
files <- files[grep("all",files)]
files <- files[grep("tsv",files)]


pdf(paste0("Plots/5_Predictability/Age/Age_MaxSubset/GSEA_",
           type,"_slopes_0.2.pdf"))

for(t in unique(sapply(files, function(x) strsplit(x,"_")[[1]][1]))){
  
  gsea <- sapply(files[grep(t,files)],
                 function(f) read_delim(paste0("Outputs/5_Predictability/Age/Age_MaxSubset/SmoothedSlopes/GSEA/Outputs/",type,"/",f),
                                        delim = "\t", escape_double = FALSE, trim_ws = TRUE), simplify = F)
  gsea <- do.call(rbind, gsea)
  
  gsea <- subset(gsea, `FDR q-val` < 0.1)
  gsea$NES <- as.numeric(gsea$NES)
  gsea <- gsea[order(gsea$NES),]
  gsea$NAME <- factor(as.character(gsea$NAME),
                      levels = gsea$NAME)
  
  print(ggplot(gsea) +
    geom_bar(aes(x = NAME, y = NES, fill = NES), stat = "identity") +
    guides(alpha = "none") + xlab("") + ylab("Enrichment score") +
    scale_fill_viridis_c(option = "cividis", direction = -1, limits = c(-2,2)) +
    ggtitle(tissue_labels_dash[t]) + coord_flip() + theme_classic() +
    theme(text = element_text(size = 20)))
}

dev.off()




# Heatmap of enrichments

terms <- c()
for(t in unique(sapply(files, function(x) strsplit(x,"_")[[1]][1]))){
  
  g <- sapply(files[grep(t,files)],
              function(f) read_delim(paste0("Outputs/5_Predictability/Age/Age_MaxSubset/SmoothedSlopes/GSEA/Outputs/",type,"/",f),
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE), simplify = F)
  g <- do.call(rbind, g)
  
  terms <- c(terms, g$NAME[g$`FDR q-val` < 0.05])
}
rm(g); gc()


gsea_mat <- matrix(nrow = length(unique(terms)), ncol = length(MAIN_TISSUES_GSEA),
                   dimnames = list(unique(terms), names(MAIN_TISSUES_GSEA)))
gsea_pval_mat <- matrix(nrow = length(unique(terms)), ncol = length(MAIN_TISSUES_GSEA),
                        dimnames = list(unique(terms), names(MAIN_TISSUES_GSEA)))
for(t in unique(sapply(files, function(x) strsplit(x,"_")[[1]][1]))){
  
  g <- sapply(files[grep(t,files)],
              function(f) read_delim(paste0("Outputs/5_Predictability/Age/Age_MaxSubset/SmoothedSlopes/GSEA/Outputs/",type,"/",f),
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE), simplify = F)
  g <- do.call(rbind, g)
  
  gsea_mat[,
           names(MAIN_TISSUES_GSEA)[MAIN_TISSUES_GSEA == t]] <- as.numeric(g[match(rownames(gsea_mat),
                                                                                   g$NAME),]$NES)
  gsea_pval_mat[,
           names(MAIN_TISSUES_GSEA)[MAIN_TISSUES_GSEA == t]] <- as.numeric(g[match(rownames(gsea_pval_mat),
                                                                                   g$NAME),]$`FDR q-val`)
}


colnames(gsea_mat) <- tissue_labels_dash[colnames(gsea_mat)]

pdf(paste0("Plots/5_Predictability/Age/Age_MaxSubset/GSEA_",
           type,"_slopes_0.2_heatmap.pdf"), height = 12)
Heatmap(t(gsea_mat), col = colorRamp2(breaks = c(-max(abs(gsea_mat), na.rm = T),
                                                 0, max(abs(gsea_mat), na.rm = T)),
                                   hcl_palette = "cividis", reverse = T),
        heatmap_legend_param = list(title = "Enrichment score"),
        row_title = "Modules", row_names_gp = gpar(fontsize = 16),
        column_names_gp = gpar(fontsize = 16))
dev.off()


pdf("Plots/Figures/Parts/module_predictability_heatmap_alloriginalslopes.pdf",
    height = 5.5, width = 9)
Heatmap(t(gsea_mat), col = colorRamp2(breaks = c(-max(abs(gsea_mat),na.rm = T),
                                                 0,max(abs(gsea_mat),na.rm = T)),
                                   # hcl_palette = "cividis", reverse = T),
                                   colors = c(age_palette["60-69"],"white",age_palette["20-29"])),
        heatmap_legend_param = list(title = "Enrichment score"),
        row_names_gp = gpar(fontsize = 16),
        column_names_gp = gpar(fontsize = 16), row_dend_side = "right",
        row_names_side = "left", column_names_rot = 40)
dev.off()

stars.pval(gsea_pval_mat)


# GO plot

pdf("Plots/Figures/Parts/module_predictability_heatmap_alloriginalslopes_GO.pdf",
    height = 11, width = 8)
Heatmap(gsea_mat, col = colorRamp2(breaks = c(-max(abs(gsea_mat),na.rm = T),
                                              0,max(abs(gsea_mat),na.rm = T)),
                                      # hcl_palette = "cividis", reverse = T),
                                      colors = c(age_palette["60-69"],"white",age_palette["20-29"])),
        heatmap_legend_param = list(title = "Enrichment score"),
        row_names_gp = gpar(fontsize = 16),
        column_names_gp = gpar(fontsize = 16), row_dend_side = "right",
        row_names_side = "left", column_names_rot = 40)
dev.off()


# Plot specific modules

net <- ReadRDS("Outputs/0_Preprocessing/stabsel_filtered_trans_largestCC_network_Hs.rds")
net <- MatrixToSquare(net)
net_graph <- AdjToIgraph(net)

# membership <- ReadRDS("Outputs/2_Clustering/stabsel_filtered_trans_largestCC_Adjacency_undirected_weightssum_rownormalized_reclustered_membership.rds")
# 
# MODULE = "6"


for(tissue in MAIN_TISSUES){
  
  slopes <- ReadRDS(paste0("Outputs/5_Predictability/Age/Age_MaxSubset/",
                           tissue,"_ageslope_well_predicted.rds"))
  s <- slopes[order(slopes[,"pval"]),]
  s <- s[1:100,]
  s <- s[order(abs(s[,"Slope"]), decreasing = T),]
  
  g <- rownames(s)[1:20]
  nbs <- ego(net_graph, order = 1, nodes = g, mode = "all", mindist = 0)
  current_graph <- induced_subgraph(net_graph,
                                    # vids = which(V(net_graph)$name %in% g))
                                    vids = unlist(nbs))
  # current_graph <- set_vertex_attr(current_graph, name = "Slope",
  #                                  value = ?)
  # V(current_graph)$color <- color_legend[V(current_graph)$Submodule]
  E(current_graph)$width <- abs(E(current_graph)$weight)*10
  g <- intersect(rownames(slopes),V(current_graph)$name)
  V(current_graph)$Slope <- NA
  V(current_graph)$Slope[match(g,V(current_graph)$name)] <- slopes[g,"Slope"]
  V(current_graph)$color <-   colorRamp2(breaks = c(-max(abs(V(current_graph)$Slope), na.rm=T),
                                                    0,max(abs(V(current_graph)$Slope), na.rm=T)),
                                         colors = c(age_palette["60-69"], "white",
                                                    age_palette["20-29"]))(V(current_graph)$Slope)
  c <- components(current_graph)
  
  pdf(paste0("Plots/5_Predictability/Age/Age_MaxSubset/smoothed_slopes_",
             tissue,"_tophits_updatedcolors.pdf"),
      # height = 13, width = 13)
      height = 6, width = 6)
  for(component in 1:c$no){
    print(plot(induced_subgraph(current_graph,
                                vids = names(which(c$membership == component))),
               layout = layout_with_lgl,
               # vertex.size = 10, vertex.label.cex = 1.5,
               vertex.size = 3,
               # vertex.size = 8,
               vertex.label.cex = 1.2,
               edge.curved = .1, vertex.label.dist = 1, edge.arrow.size = .4,
               vertex.frame.color = "black", vertex.label.color = "black",
               vertex.label.family = "Helvetica",
               main = tissue_labels_dash[tissue]))
  }
  dev.off()
  
}


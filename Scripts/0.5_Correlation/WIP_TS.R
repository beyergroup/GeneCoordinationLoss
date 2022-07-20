library(reshape2)
library(ggplot2)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)
source("Scripts/functions.R")

core.complex.list <- ReadRDS("Outputs/0.5_Correlation/coreComplexes.rds")

# single cell
sc_files <- list.files("Outputs/3_TSDataPrep/Normalized/Subset",
                       pattern = "quantile")
files_SS2 <- sc_files[grep("smartseq2",sc_files)]
files_10X <- sc_files[grep("10X",sc_files)]

# ------------------------------------------------------------------------------

# Get DESeq2-normalized Tabula Sapiens for expression and variance computation

abs_expression <- list("Smart-seq2" = list(), "10X" = list())
variance <- list("Smart-seq2" = list(), "10X" = list())
within_cor <- list("Smart-seq2" = list(), "10X" = list())

for(file in files_SS2){
  
  cell_type <- gsub("alldetected_genes_","",file)
  cell_type <- gsub("_smartseq2_quantile.rds","",cell_type)
  
  data <- ReadRDS(paste0("Outputs/3_TSDataPrep/Normalized/Subset/",file))
  
  abs_expression[["Smart-seq2"]][[cell_type]] <- rowMeans(data)
  variance[["Smart-seq2"]][[cell_type]] <- apply(data, 1, var, na.rm = T)
  
  # within_cor[["Smart-seq2"]][[cell_type]] <- lapply(core.complex.list,
  #                                                   function(genes){
  #                                                     d <- data[intersect(genes,rownames(data)),]
  #                                                     d[d == 0] <- NA
  #                                                     return(cor(t(d), use = "pairwise.complete.obs"))})
  within_cor[["Smart-seq2"]][[cell_type]] <- lapply(core.complex.list,
                                                    function(genes){
                                                      d <- data[intersect(genes,rownames(data)),,drop = F]
                                                      if(sum(rowSums(d != 0) > 0) < 2){
                                                        return(NA)
                                                      } else{
                                                        d <- d[rowSums(d != 0) > 0,]
                                                        return(cor(t(d), use = "pairwise.complete.obs"))
                                                      }})
  within_cor[["Smart-seq2"]][[cell_type]] <- 
    within_cor[["Smart-seq2"]][[cell_type]][unlist(lapply(within_cor[["Smart-seq2"]][[cell_type]],
                                                          length)) > 1]
}

for(file in files_10X){
  
  cell_type <- gsub("alldetected_genes_","",file)
  cell_type <- gsub("_10X_quantile.rds","",cell_type)
  
  data <- ReadRDS(paste0("Outputs/3_TSDataPrep/Normalized/Subset/",file))
  
  abs_expression[["10X"]][[cell_type]] <- rowMeans(data)
  variance[["10X"]][[cell_type]] <- apply(data, 1, var, na.rm = T)
  
  within_cor[["10X"]][[cell_type]] <- lapply(core.complex.list,
                                             function(genes){
                                               d <- data[intersect(genes,rownames(data)),,drop = F]
                                               if(sum(rowSums(d != 0) > 0) < 2){
                                                 return(NA)
                                               } else{
                                                 d <- d[rowSums(d != 0) > 0,]
                                                 return(cor(t(d), use = "pairwise.complete.obs"))
                                               }
                                               })
  within_cor[["10X"]][[cell_type]] <- 
    within_cor[["10X"]][[cell_type]][unlist(lapply(within_cor[["10X"]][[cell_type]],
                                                   length)) > 1]
}

saveRDS(abs_expression, "Outputs/0.5_Correlation/avg_absexpression_TS.rds")
saveRDS(variance, "Outputs/0.5_Correlation/expressionvar_TS.rds")
saveRDS(within_cor, "Outputs/0.5_Correlation/withincor_corecomplexes_TS.rds")

# ------------------------------------------------------------------------------

# Overall correlation dist of core complexes

plot.data <- lapply(within_cor,
                    function(l) lapply(l,
                                       function(ll) lapply(ll, function(mat) {
                                         if(ncol(mat) > 2){
                                           return(mean(mat[upper.tri(mat)],
                                                       na.rm = TRUE))
                                         } else{
                                           return(NA)}})))

plot.data <- rbind.data.frame(cbind.data.frame(melt(plot.data$`Smart-seq2`,
                                                    value.name = "Correlation"),
                                               "Method" = "Smart-seq2"),
                              cbind.data.frame(melt(plot.data$`10X`,
                                                    value.name = "Correlation"),
                                               "Method" = "10X"))
plot.data <- plot.data[!is.na(plot.data$Correlation),]
plot.data <- droplevels.data.frame(plot.data)

cell_types <- c("Lung_type_ii_pneumocyte","Muscle_mesenchymal_stem_cell",
                "Spleen_plasma_cell")

pdf("Plots/0.5_Correlation/corecomplex_meancor_TS.pdf", width = 55, height = 25)
ggplot(plot.data) +
  geom_boxplot(aes(x = L2, y = Correlation), color = "grey",
               fill = "white", outlier.size = 0) +
  geom_point(data = subset(plot.data, L1 %in% c(cell_types,"Merged_merged")),
             aes(x = L2, y = Correlation, color = L1, shape = Method)) +
  scale_color_manual(values = c("Merged_merged" = "black", cell_type_palette),
                     name = "Cell Type", labels = cell_type_labels) +
  xlab("Core complexes") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 20), legend.position = "top")
ggplot(subset(plot.data, Method == "Smart-seq2")) +
  geom_boxplot(aes(x = L2, y = Correlation), color = "grey",
               fill = "white", outlier.size = 0) +
  geom_point(data = subset(plot.data, L1 %in% c(cell_types,"Merged_merged")),
             aes(x = L2, y = Correlation, color = L1, shape = Method)) +
  scale_color_manual(values = c("Merged_merged" = "black", cell_type_palette),
                     name = "Cell Type", labels = cell_type_labels) +
  xlab("Core complexes") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 20), legend.position = "top")
dev.off()

# ------------------------------------------------------------------------------

# Plot average complex expression

# plot.data <- data.frame()
# 
# for(c in names(core.complex.list)){
#   
#   avg <- rbind.data.frame(cbind.data.frame(melt(lapply(abs_expression[["10X"]],
#                                                        function(x) mean(x[intersect(names(x),
#                                                                                     core.complex.list[[c]])])),
#                                                 value.name = "Mean"),
#                                            Method = "10X"),
#                           cbind.data.frame(melt(lapply(abs_expression[["Smart-seq2"]],
#                                                        function(x) mean(x[intersect(names(x),
#                                                                                     core.complex.list[[c]])])),
#                                                 value.name = "Mean"),
#                                            Method = "Smart-seq2"))
#   var <- rbind.data.frame(cbind.data.frame(melt(lapply(variance[["10X"]],
#                                                        function(x) mean(x[intersect(names(x),
#                                                                                     core.complex.list[[c]])])),
#                                                 value.name = "Var"),
#                                            Method = "10X"),
#                           cbind.data.frame(melt(lapply(variance[["Smart-seq2"]],
#                                                        function(x) mean(x[intersect(names(x),
#                                                                                     core.complex.list[[c]])])),
#                                                 value.name = "Var"),
#                                            Method = "Smart-seq2"))
#   plot.data <- rbind.data.frame(plot.data,
#                                 cbind.data.frame(merge(avg,var,"L1"),
#                                                  "Complex" = c))
# }
# 
# ggplot(plot.data) +
#   geom_density(aes(x = Mean, color = Method.x)) +
#   facet_wrap(~ L1)
# ggplot(plot.data) +
#   geom_density(aes(x = Var, color = Method.y)) +
#   facet_wrap(~ L1)
# # bc these are core complexes, avg expression is mid-high and var is low
# # Smart-seq2 data looks better
# 
# pdf("Plots/0.5_Correlation/corecomplex_meanexpr_TS.pdf", width = 55, height = 25)
# ggplot(plot.data) +
#   geom_point(data = plot.data,
#              aes(x = Complex, y = Mean, color = L1, shape = Method.x)) +
#   scale_color_manual(values = cell_type_palette, name = "Cell Type") +
#   xlab("Core complexes") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#         text = element_text(size = 20), legend.position = "top")
# dev.off()

# ------------------------------------------------------------------------------

complex = "Ribosome, cytoplasmic"
tissue = "Merged_merged"

for(method in c("Smart-seq2","10X")){
  
  h <- Heatmap(within_cor[[method]][[tissue]][[complex]], col = inferno(20),
               # row_title = paste0("Genes in complex\n(",complex,")"),
               # column_title = paste0("Genes in complex\n(",complex,")"),
               row_title = "Ribosomal genes", column_title = "Ribosomal genes",
               show_row_names = FALSE, show_column_names = FALSE,
               show_row_dend = FALSE, show_column_dend = FALSE,
               border_gp = gpar(col = "black"),
               heatmap_legend_param = list(title = "Pearson's r",
                                           direction = "horizontal",
                                           title_position = "topcenter",
                                           legend_width = unit(4,"cm")))
  heat.plot <- grid.grabExpr(draw(h, heatmap_legend_side = "bottom"))
  saveRDS(heat.plot, paste0("Outputs/0.5_Correlation/heatmap_",tissue,
                            "_",method,"_TS_ribosome_plot.rds"))
}


genes <- c("RPLP0","RPLP1")

plot.data <- data.frame()

for(file in files_SS2){
  
  cell_type <- gsub("alldetected_genes_","",file)
  cell_type <- gsub("_smartseq2_quantile.rds","",cell_type)
  
  file <- gsub("_quantile","_centered",file)
  file <- gsub("alldetected","quantified",file)
  
  data <- readRDS(paste0("Outputs/3_TSDataPrep/Centered/",
                         file))
  
  plot.data <- rbind.data.frame(plot.data, data.frame(t(as.matrix(data[genes,])),
                                                      "CellType" = cell_type,
                                                      "Method" = "Smart-seq2"))
}
for(file in files_10X){
  
  cell_type <- gsub("alldetected_genes_","",file)
  cell_type <- gsub("_10X_quantile.rds","",cell_type)
  
  file <- gsub("_quantile","_centered",file)
  file <- gsub("alldetected","quantified",file)
  
  data <- readRDS(paste0("Outputs/3_TSDataPrep/Centered/",
                         file))
  
  plot.data <- rbind.data.frame(plot.data, data.frame(t(as.matrix(data[genes,])),
                                                      "CellType" = cell_type,
                                                      "Method" = "10X"))
}

ggplot(subset(plot.data, Method == "Smart-seq2")) +
  geom_point(aes_string(x = genes[1], y = genes[2],
                        color = "CellType"))


# ------------------------------------------------------------------------------

for(method in c("10X","Smart-seq2")){
  lm.plot <- ggplot() +
    geom_smooth(data = subset(plot.data,
                              (!(CellType %in% c(cell_types,
                                                 "Merged_merged"))) &
                                (Method == method)),
                aes_string(x = genes[1], y = genes[2], group = "CellType"),
                method = "lm", se = FALSE, color = "grey") +
    geom_smooth(data = subset(plot.data,
                              (CellType %in% c(cell_types,
                                               "Merged_merged")) &
                                (Method == method)),
                aes_string(x = genes[1], y = genes[2], color = "CellType"),
                method = "lm", size = 2) +
    scale_color_manual(values = c("Merged_merged" = "black", cell_type_palette),
                       labels = cell_type_labels,
                       name = NULL) + ggtitle(paste(method,"protocol")) +
    xlab(paste(genes[1],"centered expr.")) +
    ylab(paste(genes[2],"centered expr.")) +
    theme_classic() + theme(text = element_text(size = 20),
                            legend.position = "bottom",
                            legend.text = element_text(size = 13))
  saveRDS(lm.plot, paste0("Outputs/0.5_Correlation/trendlines_",method,"_TS_",
                          genes[1],"_",genes[2],"_ribosome_plot.rds"))
}

# libs and all that crap
library(readr)
library(dorothea)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ComplexHeatmap)
library(viridis)
source("Scripts/functions.R")
source("Scripts/0.5_Correlation/params.R")

# 1. Define gene groups --------------------------------------------------------

# CORUM core complexes
coreComplexes <- read_delim("Data/CORUM/coreComplexes.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
nrow(coreComplexes) # 3512

# Subset to human
coreComplexes <- subset(coreComplexes, Organism == "Human")
nrow(coreComplexes) # 2417

# Extract gene names per complex
core.complex.list <- c()
for(i in 1:nrow(coreComplexes)){
  core.complex.list[[coreComplexes$ComplexName[i]]] <- strsplit(coreComplexes$`subunits(Gene name)`[i],
                                                                split = ";")[[1]]
}
core.complex.list <- core.complex.list[lapply(core.complex.list, length) > 5] # 346

# dorothea

tf_targets <- dorothea::dorothea_hs
tf_targets <- subset(tf_targets, confidence %in% c("A","B"))
tf.list <- sapply(unique(tf_targets$tf),
                  function(x) subset(tf_targets, tf == x)$target)

gene.sets <- c(core.complex.list,tf.list)
saveRDS(gene.sets, "Outputs/0.5_Correlation/gene_sets.rds")

rm(coreComplexes,core.complex.list,tf_targets,tf.list,i); gc()


# 2. Plot average mean expression and variance ---------------------------------

# In GTEx

gtex_expression <- readRDS("Outputs/0.5_Correlation/avg_absexpression.rds")
gtex_variance <- readRDS("Outputs/0.5_Correlation/expressionvar.rds")

plot.data <- data.frame()

for(c in names(gene.sets)){
  
  avg <- melt(lapply(gtex_expression,
                     function(x) mean(x[intersect(names(x),
                                                  gene.sets[[c]])])),
              value.name = "Mean")
  var <- melt(lapply(gtex_variance,
                     function(x) mean(x[intersect(names(x),
                                                  gene.sets[[c]])])),
              value.name = "Var")
  cor <- melt(lapply(within_cor[["GTEx"]], function(l){ 
    if(ncol(l[[c]]) >= 5){
      return(mean(l[[c]][upper.tri(l[[c]])]))
    } else{
      return(NA)}}), value.name = "Correlation")
  
  plot.data <- rbind.data.frame(plot.data,
                                cbind.data.frame(merge(cor,merge(avg,var,"L1"),
                                                       "L1"),
                                                 "Gene set" = c,
                                                 "Dataset" = "GTEx"))
}

# In Tabula Sapiens

TS_expression <- readRDS("Outputs/0.5_Correlation/avg_absexpression_TS.rds")
TS_variance <- readRDS("Outputs/0.5_Correlation/expressionvar_TS.rds")

for(method in c("10X","Smart-seq2")){
  for(c in names(gene.sets)){
    
    avg <- melt(lapply(TS_expression[[method]],
                       function(x) mean(x[intersect(names(x),
                                                    gene.sets[[c]])])),
                value.name = "Mean")
    var <- melt(lapply(TS_variance[[method]],
                       function(x) mean(x[intersect(names(x),
                                                    gene.sets[[c]])])),
                value.name = "Var")
    cor <- melt(lapply(within_cor[[method]], function(l){ 
      if(is.null(l[[c]])){
        return(NA)
      } else if(ncol(l[[c]]) >= 5){
        return(mean(l[[c]][upper.tri(l[[c]])]))
      } else{
        return(NA)}}), value.name = "Correlation")
    
    plot.data <- rbind.data.frame(plot.data,
                                  cbind.data.frame(merge(cor,merge(avg,var,
                                                                   "L1"),"L1"),
                                                   "Gene set" = c,
                                                   "Dataset" = method))
  }
}

plot.data$Dataset <- factor(as.character(plot.data$Dataset),
                            levels = c("GTEx","10X","Smart-seq2"))
plot.data$L1[plot.data$L1 %in% c("CrossTissue","Merged_merged")] <- 
  "Merged Tissues / Cell types"
plot.data$L1 <- factor(as.character(plot.data$L1),
                       levels = c("Uterus_fibroblast",
                                  "Tongue_basal_cell","Thyroid",
                                  "Spleen_plasma_cell",
                                  "Skin-SunExposed(Lowerleg)",
                                  "Skin-NotSunExposed(Suprapubic)",
                                  "Prostate_epithelial_cell",
                                  "Pancreas_pancreatic_acinar_cell",
                                  "Nerve-Tibial","Muscle_mesenchymal_stem_cell",
                                  "Muscle_endothelial_cell_of_vascular_tree",
                                  "Muscle-Skeletal",
                                  "Mammary_luminal_epithelial_cell_of_mammary_gland",
                                  "Lung_type_ii_pneumocyte","Lung_basal_cell",
                                  "Lung","Esophagus-Mucosa","Brain",
                                  "Bladder_fibroblast",
                                  "Bladder_bladder_urothelial_cell",
                                  "WholeBlood","Artery-Tibial",
                                  "Adipose-Subcutaneous",
                                  "Merged Tissues / Cell types"))
# Create Figure S1:

mean.plot <- ggplot(plot.data) +
  geom_violin(aes(y = L1, x = Mean), fill = "grey", draw_quantiles = 0.5) +
  geom_point(data = subset(plot.data, `Gene set` %in% gene_sets),
             aes(y = L1, x = Mean, color = `Gene set`), size = 2) +
  facet_wrap(~ Dataset, scales = "free_x") +
  xlab("Average expression of gene set") + ylab("Tissue / Cell type") +
  scale_y_discrete(labels = c(tissue_labels_dash,cell_type_labels_long)) +
  scale_color_manual(values = gene_set_colors[gene_sets],
                     labels = gene_set_labels[gene_sets]) +
  theme(text = element_text(size = 20), legend.position = "bottom")

var.plot <- ggplot(plot.data) +
  geom_violin(aes(y = L1, x = Var), fill = "grey", draw_quantiles = 0.5) +
  geom_point(data = subset(plot.data, `Gene set` %in% gene_sets),
             aes(y = L1, x = Var, color = `Gene set`), size = 2) +
  facet_wrap(~ Dataset, scales = "free_x") +
  xlab("Average variance of gene set") + ylab("Tissue / Cell type") +
  scale_y_discrete(labels = c(tissue_labels_dash,cell_type_labels_long)) +
  scale_color_manual(values = gene_set_colors[gene_sets],
                     labels = gene_set_labels[gene_sets]) +
  theme(text = element_text(size = 20), legend.position = "bottom")

pdf("Plots/Figures/FigureS1.pdf", width = 12, height = 20)
ggarrange(mean.plot, var.plot,
          nrow = 2, labels = "AUTO",
          align = "hv", common.legend = T,
          legend = "bottom", font.label = list(size = 22))
dev.off()

# Create panel A of Figure 1:

mean.dens.plot <- ggplot() +
  geom_density(data = subset(plot.data, (L1 == "Merged Tissues / Cell types") &
                               (Dataset %in% c("GTEx","Smart-seq2"))),
               aes(x = Mean, y = ..scaled..), outline.type = "both") +
  geom_segment(data = subset(plot.data, (L1 == "Merged Tissues / Cell types") &
                             (`Gene set` %in% gene_sets) &
                             (Dataset %in% c("GTEx","Smart-seq2"))),
             aes(x = Mean, xend = Mean, color = `Gene set`),
             y = -0.25, yend = 0.25) +
  geom_label_repel(data = subset(plot.data, (L1 == "Merged Tissues / Cell types") &
                                  (`Gene set` %in% gene_sets) &
                                  (Dataset %in% c("GTEx","Smart-seq2"))),
                  aes(x = Mean, color = `Gene set`, label = `Gene set`),
                  y = 0.35, direction = "x", fontface = 2, seed = 6,
                  fill = "white", label.size = 0, label.padding = .1) +
  facet_grid(Dataset ~ ., switch = "y", scales = "free",
             labeller = labeller(Dataset = c("GTEx" = "GTEx",
                                             "Smart-seq2" = "Tabula Sapiens"))) + 
  xlab("Average expression of gene set") + theme_classic() +
  guides(color = "none") +
  scale_color_manual(values = gene_set_colors, labels = gene_set_labels) +
  theme(text = element_text(size = 20), legend.position = "bottom",
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        strip.background = element_rect(color = "transparent",
                                        fill = "white"))

var.dens.plot <- ggplot() +
  geom_density(data = subset(plot.data, (L1 == "Merged Tissues / Cell types") &
                               (Dataset %in% c("GTEx","Smart-seq2"))),
               aes(x = Var, y = ..scaled..), outline.type = "both") +
  geom_segment(data = subset(plot.data, (L1 == "Merged Tissues / Cell types") &
                               (`Gene set` %in% gene_sets) &
                               (Dataset %in% c("GTEx","Smart-seq2"))),
               aes(x = Var, xend = Var, color = `Gene set`),
               y = -0.25, yend = 0.25) +
  geom_label_repel(data = subset(plot.data, (L1 == "Merged Tissues / Cell types") &
                                   (`Gene set` %in% gene_sets) &
                                   (Dataset %in% c("GTEx","Smart-seq2"))),
                   aes(x = Var, color = `Gene set`, label = `Gene set`),
                   y = 0.35, direction = "x", fontface = 2, seed = 5,
                   fill = "white", label.size = 0, label.padding = .1) +
  facet_grid(Dataset ~ ., scales = "free", switch = "y",
             labeller = labeller(Dataset = c("GTEx" = "GTEx",
                                             "Smart-seq2" = "Tabula Sapiens"))) + 
  xlab("Average variance of gene set") + theme_classic() +
  guides(color = "none") +
  scale_color_manual(values = gene_set_colors, labels = gene_set_labels) +
  theme(text = element_text(size = 20), legend.position = "bottom",
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        strip.background = element_rect(color = "transparent",
                                        fill = "white"))

cor.dens.plot <- ggplot() +
  geom_density(data = subset(plot.data, (L1 == "Merged Tissues / Cell types") &
                               (Dataset %in% c("GTEx","Smart-seq2"))),
               aes(x = Correlation, y = ..scaled..), outline.type = "both") +
  geom_segment(data = subset(plot.data, (L1 == "Merged Tissues / Cell types") &
                               (`Gene set` %in% gene_sets) &
                               (Dataset %in% c("GTEx","Smart-seq2"))),
               aes(x = Correlation, xend = Correlation, color = `Gene set`),
               y = -0.25, yend = 0.25) +
  geom_label_repel(data = subset(plot.data, (L1 == "Merged Tissues / Cell types") &
                                   (`Gene set` %in% gene_sets) &
                                   (Dataset %in% c("GTEx","Smart-seq2"))),
                   aes(x = Correlation, color = `Gene set`, label = `Gene set`),
                   y = 0.35, direction = "x", fontface = 2, seed = 5,
                   fill = "white", label.size = 0, label.padding = .1) +
  facet_grid(Dataset ~ ., scales = "free", switch = "y",
             labeller = labeller(Dataset = c("GTEx" = "GTEx",
                                             "Smart-seq2" = "Tabula Sapiens"))) + 
  xlab("Average correlation within gene set") + theme_classic() +
  guides(color = "none") + scale_x_continuous(limits = c(-1,1)) +
  scale_color_manual(values = gene_set_colors, labels = gene_set_labels) +
  theme(text = element_text(size = 20), legend.position = "bottom",
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        strip.background = element_rect(color = "transparent",
                                        fill = "white"))

saveRDS(list("mean" = mean.dens.plot, "var" = var.dens.plot),
        "Outputs/0.5_Correlation/density_plots_crosstissue.rds")



# 3. Within-set correlation ----------------------------------------------------

within_cor <- list("GTEx" = list(), "10X" = list(), "Smart-seq2" = list())

# GTEx

gtex_files <- list.files("GTEx_Networks/Tissue_Networks/Outputs",
                         pattern = "sampled_data.rds", full.names = T)

for(file in gtex_files){
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  data <- readRDS(file)
  data <- DataENSGToSymbol(data, remove_dup = T)
  
  within_cor[["GTEx"]][[tissue]] <-
    lapply(gene.sets, function(genes) cor(t(data[intersect(genes,
                                                           rownames(data)),])))
  within_cor[["GTEx"]][[tissue]][["Random"]] <- cor(t(data[intersect(random_genes,
                                                                     rownames(data)),]))
}

# Tabula Sapiens

sc_files <- list.files("Outputs/3_TSDataPrep/Normalized/Subset",
                       pattern = "quantile")
files_SS2 <- sc_files[grep("smartseq2",sc_files)]
files_10X <- sc_files[grep("10X",sc_files)]
rm(sc_files); gc()

for(file in files_SS2){
  
  cell_type <- gsub("alldetected_genes_","",file)
  cell_type <- gsub("_smartseq2_quantile.rds","",cell_type)
  
  data <- ReadRDS(paste0("Outputs/3_TSDataPrep/Normalized/Subset/",file))
  
  within_cor[["Smart-seq2"]][[cell_type]] <- lapply(gene.sets,
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
                                                          length)) > 25]
  
  d <- data[intersect(random_genes,rownames(data)),,drop = F]
  if(sum(rowSums(d != 0) > 0) < 2){
    within_cor[["Smart-seq2"]][[cell_type]][["Random"]] <- NA
  } else{
    d <- d[rowSums(d != 0) > 0,]
    within_cor[["Smart-seq2"]][[cell_type]][["Random"]] <- 
      cor(t(d), use = "pairwise.complete.obs")
  }
  rm(d)
}

for(file in files_10X){
  
  cell_type <- gsub("alldetected_genes_","",file)
  cell_type <- gsub("_10X_quantile.rds","",cell_type)
  
  data <- ReadRDS(paste0("Outputs/3_TSDataPrep/Normalized/Subset/",file))
  
  within_cor[["10X"]][[cell_type]] <- lapply(gene.sets,
                                             function(genes){
                                               d <- data[intersect(genes,rownames(data)),,drop = F]
                                               if(sum(rowSums(d != 0) > 0) < 2){
                                                 return(NA)
                                               } else{
                                                 d <- d[rowSums(d != 0) > 0,]
                                                 return(cor(t(d), use = "pairwise.complete.obs"))
                                               }})
  within_cor[["10X"]][[cell_type]] <- 
    within_cor[["10X"]][[cell_type]][unlist(lapply(within_cor[["10X"]][[cell_type]],
                                                          length)) > 25]
  
  d <- data[intersect(random_genes,rownames(data)),,drop = F]
  if(sum(rowSums(d != 0) > 0) < 2){
    within_cor[["10X"]][[cell_type]][["Random"]] <- NA
  } else{
    d <- d[rowSums(d != 0) > 0,]
    within_cor[["10X"]][[cell_type]][["Random"]] <- 
      cor(t(d), use = "pairwise.complete.obs")
  }
  rm(d)
}
rm(files_SS2,files_10X); gc()

WriteRDS(within_cor, "Outputs/0.5_Correlation/withincor_genesets.rds")


# 4. Average within-set correlation --------------------------------------------

cor.plot.data <- data.frame()
for(method in c("Smart-seq2","GTEx")){
  cor.plot.data <- rbind.data.frame(cor.plot.data,
                                    data.frame(melt(lapply(within_cor[[method]][[grep("CrossTissue|Merged_merged",
                                                                                      names(within_cor[[method]]))]],
                                                           function(mat) {
                                                             if(ncol(mat) >= 5){
                                                               return(mean(mat[upper.tri(mat)]))
                                                             } else{
                                                               return(NA)}}), value.name = "Correlation"),
                                               "Dataset" = method))
}



ggplot(plot.data) +
  geom_density(aes(x = Correlation, y = ..scaled..), outline.type = "both") +
  geom_segment(data = subset(plot.data, L1 %in% c(gene_sets,"Random")),
               aes(x = Correlation, xend = Correlation, color = L1),
               y = -0.25, yend = 0.25) +
  geom_label_repel(data = subset(plot.data, L1 %in% c(gene_sets,"Random")),
                   aes(x = Correlation, color = L1, label = L1),
                   y = 0.35, direction = "x", fontface = 2, seed = 5,
                   fill = "white", label.size = 0, label.padding = .1) +
  facet_grid(Dataset ~ ., scales = "free", switch = "y",
             labeller = labeller(Dataset = c("GTEx" = "GTEx",
                                             "Smart-seq2" = "Tabula Sapiens"))) + 
  xlab("Average correlation within gene set") + theme_classic() +
  guides(color = "none") +
  scale_color_manual(values = gene_set_colors, labels = gene_set_labels) +
  theme(text = element_text(size = 20), legend.position = "bottom",
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        strip.background = element_rect(color = "transparent",
                                        fill = "white"))


# 5. Expression of gene set members --------------------------------------------

plot.data <- data.frame()

for(gene_set in gene_sets){
  plot.data <- rbind.data.frame(plot.data,
                                data.frame(melt(do.call(cbind,
                                                        lapply(gtex_expression,
                                                               function(m) m[gene.sets[[gene_set]]])),
                                                varnames = c("Gene","Tissue"),
                                                value.name = "Mean"),
                                           "Dataset" = "GTEx",
                                           "GeneSet" = gene_set))
  plot.data <- rbind.data.frame(plot.data,
                                data.frame(melt(do.call(cbind,
                                                        lapply(TS_expression$`Smart-seq2`,
                                                               function(m) m[gene.sets[[gene_set]]])),
                                                varnames = c("Gene","Tissue"),
                                                value.name = "Mean"),
                                           "Dataset" = "Smart-seq2",
                                           "GeneSet" = gene_set))
}

# Add random set
set.seed(1)
random_genes <- sample(setdiff(intersect(names(gtex_expression$CrossTissue),
                                         names(which(table(unlist(lapply(TS_expression$`Smart-seq2`,
                                                                         names))) == 13))),
                               unlist(gene.sets[gene_sets])), 85)
plot.data <- rbind.data.frame(plot.data,
                              data.frame(melt(do.call(cbind,
                                                      lapply(gtex_expression,
                                                             function(m) m[random_genes])),
                                              varnames = c("Gene","Tissue"),
                                              value.name = "Mean"),
                                         "Dataset" = "GTEx",
                                         "GeneSet" = "Random"))
plot.data <- rbind.data.frame(plot.data,
                              data.frame(melt(do.call(cbind,
                                                      lapply(TS_expression$`Smart-seq2`,
                                                             function(m) m[random_genes])),
                                              varnames = c("Gene","Tissue"),
                                              value.name = "Mean"),
                                         "Dataset" = "Smart-seq2",
                                         "GeneSet" = "Random"))

plot.data$GeneSet <- factor(as.character(plot.data$GeneSet),
                            levels = c(gene_sets,"Random"))
plot.data$Tissue <- factor(as.character(plot.data$Tissue),
                           levels = c(names(tissue_labels_long),
                                      names(cell_type_labels_long)))

# Create panel B of Figure 1:

expr.boxplot <- ggplot(plot.data) +
  geom_boxplot(aes(x = Tissue, y = Mean, fill = GeneSet)) +
  facet_grid(~ Dataset, scales = "free", drop = T, switch = "x",
             labeller = labeller(Dataset = c("GTEx" = "GTEx",
                                             "Smart-seq2" = "Tabula Sapiens (Smart-seq2)"))) +
  scale_fill_manual(values = gene_set_colors, labels = gene_set_labels) +
  ylab("Average expression of genes in set") + xlab("Tissues / Cell Types") +
  theme_classic() + labs(fill = guides(title = "Gene set")) +
  scale_x_discrete(labels = c(tissue_labels_long, cell_type_labels_long)) +
  theme(text = element_text(size = 20), legend.position = "bottom",
        strip.background = element_rect(color = "transparent", fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

saveRDS(expr.boxplot, "Outputs/0.5_Correlation/expr_boxplot.rds")


# 6. Correlation heatmap -------------------------------------------------------

complex = gene_sets[2]

for(method in c("GTEx","Smart-seq2")){
  
  h <- Heatmap(within_cor[[method]][[grep("CrossTissue|Merged_merged",
                                          names(within_cor[[method]]))]][[complex]],
               col = inferno(20),
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



plot.data <- data.frame()
for(dataset in c("GTEx","Smart-seq2")){
  plot.data <- rbind.data.frame(plot.data,
                                data.frame(melt(lapply(within_cor[[dataset]],
                                                       function(l) lapply(l[c(gene_sets,"Random")],
                                                                          function(ll) ll[upper.tri(ll)])),
                                                value.name = "Correlation"),
                                           "Dataset" = dataset))
}

plot.data$L2 <- factor(as.character(plot.data$L2),
                       levels = c(gene_sets,"Random"))
plot.data$L1 <- factor(as.character(plot.data$L1),
                       levels = c(names(tissue_labels_long),
                                  names(cell_type_labels_long)))


# Create panel D of Figure 1:

ggplot(plot.data) +
  geom_boxplot(aes(x = L1, y = Correlation, fill = L2)) +
  facet_grid(~ Dataset, scales = "free", drop = T, switch = "x",
             labeller = labeller(Dataset = c("GTEx" = "GTEx",
                                             "Smart-seq2" = "Tabula Sapiens (Smart-seq2)"))) +
  scale_fill_manual(values = gene_set_colors, labels = gene_set_labels) +
  ylab("Average expression of genes in set") + xlab("Tissues / Cell Types") +
  theme_classic() + labs(fill = guides(title = "Gene set")) +
  scale_x_discrete(labels = c(tissue_labels_long, cell_type_labels_long)) +
  theme(text = element_text(size = 20), legend.position = "bottom",
        strip.background = element_rect(color = "transparent", fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

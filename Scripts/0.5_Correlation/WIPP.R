# libs and all that crap
library(dorothea)
library(ggpubr)
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
  plot.data <- rbind.data.frame(plot.data,
                                cbind.data.frame(merge(avg,var,"L1"),
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
    plot.data <- rbind.data.frame(plot.data,
                                  cbind.data.frame(merge(avg,var,"L1"),
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
  scale_color_manual(values = gene_set_colors, labels = gene_set_labels) +
  theme(text = element_text(size = 20), legend.position = "bottom")

var.plot <- ggplot(plot.data) +
  geom_violin(aes(y = L1, x = Var), fill = "grey", draw_quantiles = 0.5) +
  geom_point(data = subset(plot.data, `Gene set` %in% gene_sets),
             aes(y = L1, x = Var, color = `Gene set`), size = 2) +
  facet_wrap(~ Dataset, scales = "free_x") +
  xlab("Average variance of gene set") + ylab("Tissue / Cell type") +
  scale_y_discrete(labels = c(tissue_labels_dash,cell_type_labels_long)) +
  scale_color_manual(values = gene_set_colors, labels = gene_set_labels) +
  theme(text = element_text(size = 20), legend.position = "bottom")

pdf("Plots/Figures/FigureS1.pdf", width = 12, height = 20)
ggarrange(mean.plot, var.plot,
          nrow = 2, labels = "AUTO",
          align = "hv", common.legend = T,
          legend = "bottom", font.label = list(size = 22))
dev.off()

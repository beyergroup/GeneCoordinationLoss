
library(readr)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)
source("Scripts/functions.R")

# ------------------------------------------------------------------------------

# Get protein complexes

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

saveRDS(core.complex.list, "Outputs/0.5_Correlation/coreComplexes.rds")

# ------------------------------------------------------------------------------

# Get DESeq2-normalized GTEx for expression and variance computation

gtex_files <- list.files("GTEx_Networks/Tissue_Networks/Outputs",
                         pattern = "sampled_data.rds", full.names = T)

gtex_norm <- readRDS("Data/GTEx/DESeq2_normalized_gtex.rds")
gtex_norm <- DataENSGToSymbol(gtex_norm, remove_dup = T)

abs_expression <- list()
variance <- list()

for(file in gtex_files){
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  data <- readRDS(file)
  samples <- colnames(data)
  rm(data); gc()
  
  data <- gtex_norm[,samples]
  data <- log2(data+1)
  
  abs_expression[[tissue]] <- rowMeans(data)
  variance[[tissue]] <- apply(data, 1, var, na.rm = T)
}

saveRDS(abs_expression, "Outputs/0.5_Correlation/avg_absexpression.rds")
saveRDS(variance, "Outputs/0.5_Correlation/expressionvar.rds")

# ------------------------------------------------------------------------------

# Plot average complex expression

plot.data <- data.frame()

for(c in names(core.complex.list)){
  
  avg <- melt(lapply(abs_expression,
                     function(x) mean(x[intersect(names(x),
                                                  core.complex.list[[c]])])),
              value.name = "Mean")
  var <- melt(lapply(variance,
                     function(x) mean(x[intersect(names(x),
                                                  core.complex.list[[c]])])),
              value.name = "Var")
  plot.data <- rbind.data.frame(plot.data,
                                cbind.data.frame(merge(avg,var,"L1"),
                                                 "Complex" = c))
}

ggplot(plot.data) +
  geom_density(aes(x = Mean)) +
  facet_wrap(~ L1)
# bc these are core complexes, avg expression is mid-high and var is low

# ------------------------------------------------------------------------------

# Within complex correlation

gtex_files <- list.files("GTEx_Networks/Tissue_Networks/Outputs",
                         pattern = "sampled_data.rds", full.names = T)

within_cor <- list()

for(file in gtex_files){
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  data <- readRDS(file)
  data <- DataENSGToSymbol(data, remove_dup = T)
  
  within_cor[[tissue]] <- lapply(core.complex.list,
                                 function(genes) cor(t(data[intersect(genes,
                                                                      rownames(data)),])))
}

saveRDS(within_cor, "Outputs/0.5_Correlation/withincor_corecomplexes.rds")

# ------------------------------------------------------------------------------

# Overall correlation dist of core complexes

tissues = c("Adipose-Subcutaneous","Brain","Muscle-Skeletal")

plot.data <- within_cor[tissues]
plot.data <- lapply(plot.data,
                    function(l) lapply(l,
                                       function(mat) {
                                         if(ncol(mat) > 2){
                                           return(mat[upper.tri(mat)])
                                         } else{
                                           return(NA)}}))
plot.data <- melt(plot.data, value.name = "Correlation")

ggplot(plot.data) +
  geom_violin(aes(x = L2, y = Correlation)) +
  facet_wrap(~ L1, ncol = 1) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ------------------------------------------------------------------------------

# Mean correlation within core complex across tissues

plot.data <- lapply(within_cor,
                    function(l) lapply(l,
                                       function(mat) {
                                         if(ncol(mat) >2){
                                           return(mean(mat[upper.tri(mat)]))
                                         } else{
                                           return(NA)}}))
plot.data <- melt(plot.data, value.name = "Correlation")

pdf("Plots/0.5_Correlation/corecomplex_meancor.pdf", width = 55, height = 25)
ggplot(plot.data) +
  geom_boxplot(aes(x = L2, y = Correlation), color = "grey",
               fill = "white", outlier.size = 0) +
  geom_point(data = subset(plot.data, L1 %in% c(tissues,"CrossTissue")),
             aes(x = L2, y = Correlation, color = L1)) +
  scale_color_manual(values = c("CrossTissue" = "black", tissue_palette)) +
  xlab("Core complexes") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 20))
dev.off()

# ------------------------------------------------------------------------------

# Pick one complex for illustration purposes

complex = "Ribosome, cytoplasmic"
tissues = c("Adipose-Subcutaneous","Brain","Muscle-Skeletal")

plot.data <- data.frame()

for(tissue in tissues){
  
  plot.data <- rbind.data.frame(plot.data,
        data.frame("Expression" = abs_expression[[tissue]][intersect(core.complex.list[[complex]],
                                                                     names(abs_expression[[tissue]]))],
                   "Variance" = variance[[tissue]][intersect(core.complex.list[[complex]],
                                                             names(variance[[tissue]]))],
                   "Tissue" = tissue))
}

tissue_labels <- gsub("-","\n(",tissues)
tissue_labels[grepl("(",tissue_labels,
                    fixed=T)] <- paste0(tissue_labels[grepl("\n(",tissue_labels,
                                                            fixed=T)],")")
names(tissue_labels) <- tissues

expr.plot <- ggplot(plot.data) +
  geom_boxplot(aes(x = Tissue, y = Expression, fill = Tissue),
               size = 1, outlier.colour = "transparent") +
  geom_jitter(aes(x = Tissue, y = Expression), width = 0.3) +
  scale_fill_manual(values = tissue_palette) +
  scale_x_discrete(labels = tissue_labels) +
  guides(fill = "none") + ylab("Mean expression\n(log2-transformed)") +
  theme_classic() + theme(text = element_text(size = 20))

var.plot <- ggplot(plot.data) +
  geom_boxplot(aes(x = Tissue, y = Variance, fill = Tissue),
               size = 1, outlier.colour = "transparent") +
  geom_jitter(aes(x = Tissue, y = Variance), width = 0.3) +
  scale_fill_manual(values = tissue_palette) +
  scale_x_discrete(labels = tissue_labels) +
  guides(fill = "none") + ylab("Expression variance\n(log2-transformed)") +
  theme_classic() + theme(text = element_text(size = 20))

# ------------------------------------------------------------------------------

# 

genes <- c("RPLP0","RPLP1")

gtex_files <- list.files("GTEx_Networks/Tissue_Networks/Outputs",
                         pattern = "sampled_centered_data.rds", full.names = T)

plot.data <- data.frame()

# for(file in gtex_files[grep(paste(tissues,collapse="|"),gtex_files)]){
for(file in gtex_files){
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  data <- readRDS(file)
  data <- DataENSGToSymbol(data, remove_dup = T)
  
  plot.data <- rbind.data.frame(plot.data, data.frame(t(data[genes,]),
                                                      "Tissue" = tissue))
}

ggplot(plot.data) +
  geom_point(aes_string(x = genes[1], y = genes[2])) +
  facet_wrap(~ Tissue, ncol = 1)

pdf(paste0("Plots/0.5_Correlation/",genes[1],"_",genes[2],"_tissues.pdf"),
    height = 6, width = 6)
ggplot(plot.data) +
  geom_smooth(data = subset(plot.data, !(Tissue %in% c(tissues,"CrossTissue"))),
              aes_string(x = genes[1], y = genes[2], group = "Tissue"),
              method = "lm", se = FALSE, color = "grey") +
  geom_smooth(data = subset(plot.data, Tissue %in% c(tissues,"CrossTissue")),
              aes_string(x = genes[1], y = genes[2], color = "Tissue"),
              method = "lm", size = 2) +
  scale_color_manual(values = c("CrossTissue" = "black", tissue_palette),
                     labels = tissue_labels[c("CrossTissue",
                                              names(tissue_palette))],
                     name = NULL) +
  xlab(paste(genes[1],"centered expr.")) +
  ylab(paste(genes[2],"centered expr.")) +
  theme_classic() + theme(text = element_text(size = 20),
                          legend.position = "bottom",
                          legend.text = element_text(size = 13))
dev.off()

# complex = "Ribosome, cytoplasmic"
complex = "MLL−HCF complex"
tissue = "CrossTissue"

pdf(paste0("Plots/0.5_Correlation/ribosome_heatmap_",tissue,".pdf"),
    height = 5, width = 5)
h <- Heatmap(within_cor[[tissue]][[complex]], col = inferno(20),
             row_title = paste0("Genes in complex\n(",complex,")"),
             column_title = paste0("Genes in complex\n(",complex,")"),
             show_row_names = FALSE, show_column_names = FALSE,
             show_row_dend = FALSE, show_column_dend = FALSE,
             border_gp = gpar(col = "black"),
             heatmap_legend_param = list(title = "Pearson's r",
                                         direction = "horizontal",
                                         title_position = "topcenter",
                                         legend_width = unit(4,"cm")))
# draw(h, heatmap_legend_side = "bottom")
dev.off()

heat.plot <- grid.grabExpr(draw(h, heatmap_legend_side = "bottom"))
lm.plot <- ggplot(plot.data) +
  geom_smooth(data = subset(plot.data, !(Tissue %in% c(tissues,"CrossTissue"))),
              aes_string(x = genes[1], y = genes[2], group = "Tissue"),
              method = "lm", se = FALSE, color = "grey") +
  geom_smooth(data = subset(plot.data, Tissue %in% c(tissues,"CrossTissue")),
              aes_string(x = genes[1], y = genes[2], color = "Tissue"),
              method = "lm", size = 2) +
  scale_color_manual(values = c("CrossTissue" = "black", tissue_palette),
                     labels = tissue_labels[c("CrossTissue",
                                              names(tissue_palette))],
                     name = NULL) +
  xlab(paste(genes[1],"centered expr.")) +
  ylab(paste(genes[2],"centered expr.")) +
  theme_classic() + theme(text = element_text(size = 20),
                          legend.position = "bottom",
                          legend.text = element_text(size = 13))

saveRDS(heat.plot,
        paste0("Outputs/0.5_Correlation/heatmap_",
               tissue,"_GTEx_ribosome_plot.rds"))
saveRDS(lm.plot, paste0("Outputs/0.5_Correlation/trendlines_GTEx_",
                        genes[1],"_",genes[2],"_ribosome_plot.rds"))

# ------------------------------------------------------------------------------



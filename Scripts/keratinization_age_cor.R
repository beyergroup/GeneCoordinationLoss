# Look into keratins

library(biomaRt, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(reshape2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggplot2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggpubr, lib.loc = "Resources/Rlibs/R-4.0.3/")

# GO_IDs <- c("keratinization" = "GO:0031424")
GO_IDs <- c(NULL,
            "translational initiation" = "GO:0006413", # 171 genes
            # "mitochondrial respiratory chain" = "GO:0032981", # 82 genes
            # "mitochondrial respiratory chain (II-IV)" = "GO:0005749",
            # "mitochondrial respiratory chain (II-IV)" = "GO:0005750",
            # "mitochondrial respiratory chain (II-IV)" = "GO:0005751",
            NULL) # all mito: 111 genes

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","go_id"),
               filters = "go", values = GO_IDs, mart = ensembl)
genes <- subset(genes, go_id %in% GO_IDs)
genes$GOTerm <- names(GO_IDs)[match(genes$go_id,GO_IDs)]
# ribo <- genes[grepl("RPS", genes$hgnc_symbol) & !grepl("MRP",genes$hgnc_symbol),]
genes <- subset(genes, GOTerm != "structural constituent of ribosome")
# genes <- rbind.data.frame(genes, ribo)
genes <- genes$hgnc_symbol # 240 genes

k <- readRDS("~/Downloads/Keratin CROSS.rds")
genes <- k$SYMBOL # 58 genes (53 of which are part of the keratin GO term)

# Load correlations AgeTissue

files <- list.files("Outputs/Human_Network/stabsel/Predictability/AgeTissue",
  pattern = "sampled_correlations.rds", full.names = T)

correlations <- sapply(files, readRDS)
plot.data <- melt(correlations)
plot.data$Var2 <- sapply(as.character(plot.data$Var2), function(x) strsplit(tail(strsplit(x, "/")[[1]],1),
   "_sampled_correlations")[[1]][1])
plot.data$Tissue <- sapply(as.character(plot.data$Var2), function(x) strsplit(x, "_")[[1]][1])
plot.data$Tissue <- sapply(plot.data$Tissue, function(x) gsub("(","\n(",x, fixed = T))
plot.data$Age <- sapply(plot.data$Var2, function(x) strsplit(x, "_")[[1]][2])


plots <- list()
for(tissue in unique(plot.data$Tissue)){
  plots[[tissue]] <- eval(ggplot(subset(plot.data, (Tissue == tissue) & (Var1 %in% genes))) +
    geom_violin(aes(y = value, x = Age, fill = Age), size = 1, draw_quantiles = 0.5) +
    scale_fill_viridis_d(option = "magma", begin = .3, end = .7) +
    ylim(c(-1,1)) + geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Correlation coefficient") + ggtitle(tissue) + xlab("Age") +
    theme_classic() + theme(text = element_text(size = 18),
      legend.position = "bottom", legend.title = element_blank(), title = element_text(size = 18)))
}

# pdf("Plots/Human_Network/stabsel/Predictability/GOterms/keratinization_agecors_violins_FranciscoCor.pdf", height = 9, width = 16)
# pdf("Plots/Human_Network/stabsel/Predictability/GOterms/mitochondrial_agecors_violins.pdf", height = 9, width = 16)
# pdf("Plots/Human_Network/stabsel/Predictability/GOterms/ribosomal_agecors_violins.pdf", height = 9, width = 16)
pdf("Plots/Human_Network/stabsel/Predictability/GOterms/translation_agecors_violins.pdf", height = 9, width = 16)
ggarrange(plotlist = plots, align = "hv", common.legend = T, legend = "bottom")
dev.off()

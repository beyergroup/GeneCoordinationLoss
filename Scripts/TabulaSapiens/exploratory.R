# Quality control per tissue

INDIR = "../../Outputs/Tabula_Sapiens/pertissue"
METHOD = "10X"


files <- list.files(paste0(INDIR,"/",METHOD), full.names = T)

.libPaths("../../Resources/Rlibs/R-4.0.3/")
library(Seurat)
library(ggplot2)
library(formattable)

for(file in files){
  
  sce <- readRDS(file)
  
  plot.data <- as.data.frame(sce@colData)
  
  # quantify reads from mitochondrial genes
  mt_genes <- rownames(assay(sce,"raw_counts"))[grep("MT-",rownames(assay(sce,"raw_counts")))]
  mt_genes <- mt_genes[-grep("INMT-|ASMT-",mt_genes)]
  plot.data$MTpercent <- colSums(assay(sce,"raw_counts")[mt_genes,])*100/colSums(assay(sce,"raw_counts"))
  rm(mt_genes)
  
  pdf(paste0("../../Plots/Tabula_Sapiens/QC/",unique(plot.data$organ_tissue),
             c("smartseq2" = "Smart-seq2",
               "10X" = "10X")[METHOD],".pdf"))
  
  # number of cells per cell type
  tb <- as.data.frame(table(plot.data$cell_ontology_class))
  colnames(tb) <- c("Cell Type", "Numbers")
  print(formattable(tb))
  
  # library size distribution
  print(ggplot(plot.data) +
    geom_jitter(aes(x = cell_ontology_class, y = n_counts_UMIs,
                    color = cell_ontology_class), width = 0.3) +
    geom_violin(aes(x = cell_ontology_class, y = n_counts_UMIs),
                size = 1.5, fill = "transparent", draw_quantiles = 0.5) +
    scale_y_log10() + guides(color = FALSE) +
    xlab("Cell Type") + ylab("Library size (log scale)") +
    ggtitle(paste(unique(plot.data$organ_tissue),
                  c("smartseq2" = "Smart-seq2",
                    "10X" = "10X")[METHOD])) +
    theme_classic() + theme(text = element_text(size = 20),
                            axis.text.x = element_text(angle = 90,
                                                       hjust = 1,
                                                       vjust = .5)))
  
  # number of detected features 
  print(ggplot(plot.data) +
    geom_jitter(aes(x = cell_ontology_class, y = n_genes,
                    color = cell_ontology_class), width = 0.3) +
    geom_violin(aes(x = cell_ontology_class, y = n_genes),
                size = 1.5, fill = "transparent", draw_quantiles = 0.5) +
    guides(color = FALSE) +
    xlab("Cell Type") + ylab("# detected genes") +
    ggtitle(paste(unique(plot.data$organ_tissue),
                  c("smartseq2" = "Smart-seq2",
                    "10X" = "10X")[METHOD])) +
    theme_classic() + theme(text = element_text(size = 20),
                            axis.text.x = element_text(angle = 90,
                                                       hjust = 1,
                                                       vjust = .5)))
  
  # %MT features
  print(ggplot(plot.data) +
    geom_jitter(aes(x = cell_ontology_class, y = MTpercent,
                    color = cell_ontology_class), width = 0.3) +
    geom_violin(aes(x = cell_ontology_class, y = MTpercent),
                size = 1.5, fill = "transparent", draw_quantiles = 0.5) +
    guides(color = FALSE) + ylim(c(0,100)) +
    xlab("Cell Type") + ylab("% mitochondrial features") +
    ggtitle(paste(unique(plot.data$organ_tissue),
                  c("smartseq2" = "Smart-seq2",
                    "10X" = "10X")[METHOD])) +
    theme_classic() + theme(text = element_text(size = 20),
                            axis.text.x = element_text(angle = 90,
                                                       hjust = 1,
                                                       vjust = .5)))
  
  # scatter plot
  print(ggplot(plot.data) +
    geom_point(aes(y = n_genes, x = n_counts_UMIs,
                   color = MTpercent)) +
    scale_x_log10() +
    scale_color_viridis_c(name = "% MT") +
    xlab("Library size (log scale)") +
    ylab("# detected genes") +
    ggtitle(paste(unique(plot.data$organ_tissue),
                  c("smartseq2" = "Smart-seq2",
                    "10X" = "10X")[METHOD])) +
    theme_classic() + theme(text = element_text(size = 20)))
  
  dev.off()
}

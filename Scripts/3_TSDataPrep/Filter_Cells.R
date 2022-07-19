# Quality control of cells
# Apply filters on cell library sizes and mitochondrial reads percentages:
# %MT not above 4 median absolute deviations of the median across cells
# cell library sizes within an interval of 3MAD of the median across cells (in log10 space)

args <- commandArgs(trailingOnly = T)
WDIR = args[1]
METHOD = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
source("Scripts/functions.R")
source(paste0("Scripts/",WDIR,"/params.R"))
library(SingleCellExperiment)
library(ggplot2)
library(cowplot)
library(tidyverse)


files <- list.files(paste0("Outputs/3_TSDataPrep/",METHOD),
                    pattern = ".rds", full.names = T)


for(file in files){
  
  sce <- ReadRDS(file)
  
  tissue <- as.character(unique(sce$organ_tissue))
  
  # Exclude cell types with less than 100 cells
  cell_types <- names(which(table(sce$cell_ontology_class) > NCELL_THRE[METHOD]))
  if(length(cell_types) == 0)
    next
  sce <- sce[, sce$cell_ontology_class %in% cell_types]
  
  # quantify reads from mitochondrial genes
  mt_genes <- rownames(assay(sce,"raw_counts"))[grep("MT-",rownames(assay(sce,"raw_counts")))]
  mt_genes <- mt_genes[-grep("INMT-|ASMT-",mt_genes)]
  sce$MTpercent <- colSums(assay(sce,"raw_counts")[mt_genes,])*100/colSums(assay(sce,"raw_counts"))
  rm(mt_genes)
  
  pdf(paste0("Plots/",WDIR,"/Filters/cellfilters_",tissue,"_",METHOD,".pdf"),
      height = 9)
  
  # Apply filters per cell type
  
  keep <- c()
  
  for(cell_type in cell_types){
    
    message(cell_type)
    
    plots <- list()
    
    tmp_sce <- sce[, sce$cell_ontology_class == cell_type]
    
    # compute threshold
    med <- median(log10(tmp_sce$n_counts_UMIs))
    md <- mad(log10(tmp_sce$n_counts_UMIs))
    doublet_thre <- 10^(med + (3*md))
    dead_thre <- 10^(med - (3*md))
    rm(med,md); gc()
    med <- median(tmp_sce$MTpercent)
    md <- mad(tmp_sce$MTpercent)
    MT_thre <- med + (4*md)
    rm(med,md); gc()
    
    # plot
    plots[[1]] <- ggplot(as.data.frame(tmp_sce@colData)) +
      geom_point(aes(x = n_counts_UMIs, y = n_genes, color = MTpercent)) +
      geom_vline(xintercept = c(doublet_thre,dead_thre),
                 linetype = "dashed", col = "tomato") +
      scale_color_viridis_c(name = "% MT") +
      scale_x_log10() +
      xlab("Library size (log scale)") +
      ylab("# detected genes") +
      ggtitle(paste0(cell_type," - ",METHOD),tissue) +
      theme(text = element_text(size = 20))
    plots[[2]] <- ggplot(as.data.frame(tmp_sce@colData)) +
      geom_density(aes(x = n_counts_UMIs)) +
      geom_vline(xintercept = c(doublet_thre,dead_thre),
                 linetype = "dashed", col = "tomato") +
      xlab("Library size (log scale)") +
      scale_x_log10() +
      theme(text = element_text(size = 20))
    plots[[3]] <- ggplot(as.data.frame(tmp_sce@colData)) +
      geom_density(aes(x = MTpercent)) +
      geom_vline(xintercept = MT_thre,
                 linetype = "dashed", col = "violetred") +
      xlab("% Mitochondrial features") +
      theme(text = element_text(size = 20))
    print(plot_grid(plotlist = plots, ncol = 1, align = "v", axis = "lr",
                    rel_heights = c(2,1,1)))
    
    # apply to data
    keep <- c(keep, names(which((tmp_sce$n_counts_UMIs <= doublet_thre) &
                                  (tmp_sce$MTpercent <= MT_thre) &
                                  (tmp_sce$n_counts_UMIs) >= dead_thre)))
    
    rm(tmp_sce,plots,doublet_thre,MT_thre,dead_thre); gc()
  }
  
  dev.off()
  
  sce <- sce[, keep]
  WriteRDS(sce, paste0("Outputs/",WDIR,"/Filters/cellfiltered_",
                      METHOD,"_",tissue,".rds"))
  
  rm(sce,tissue,cell_types,cell_type,keep); gc()
}

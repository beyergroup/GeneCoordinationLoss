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

GSEA_DIR = "Outputs/5_Predictability/WellPredicted_TissueFilters/GSEA/Out"
TYPE = "GOmf"
PLOT_DIR = "Plots/5_Predictability/WellPredicted_TissueFilters"

files <- list.files(paste0(GSEA_DIR,"/",TYPE),
                    recursive = T, pattern = "gsea_report_for_na_")
files <- files[grep("tsv",files)]


pdf(paste0(PLOT_DIR,"/GSEA_",
           TYPE,"_predictability.pdf"))

for(t in unique(sapply(files, function(x) strsplit(x,"\\.")[[1]][1]))){
  
  gsea <- sapply(files[grep(t,files)],
                 function(f) read_delim(paste0(GSEA_DIR,"/",TYPE,"/",f),
                                        delim = "\t", escape_double = FALSE, trim_ws = TRUE), simplify = F)
  gsea <- do.call(rbind, gsea)
  
  gsea$NES <- as.numeric(gsea$NES)
  gsea <- subset(gsea, (`FDR q-val` < 0.1) & (abs(NES) > 1.5))
  gsea <- gsea[order(gsea$NES),]
  gsea$NAME <- factor(as.character(gsea$NAME),
                      levels = gsea$NAME)
  
  print(ggplot(gsea) +
    geom_bar(aes(x = NAME, y = NES, fill = NES), stat = "identity") +
    guides(alpha = "none") + xlab("") + ylab("Enrichment score") +
    scale_fill_gradient2(low = "#FF6347", mid = "white", high = "#40547c") +
    ggtitle(tissue_labels_dash_all[names(which(ALL_TISSUES_GSEA == t))]) +
    coord_flip() + theme_classic() +
    theme(text = element_text(size = 20)))
}

dev.off()




# Heatmap of enrichments

terms <- c()
for(t in unique(sapply(files, function(x) strsplit(x,"\\.")[[1]][1]))){
  
  g <- sapply(files[grep(t,files)],
              function(f) read_delim(paste0(GSEA_DIR,"/",TYPE,"/",f),
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE), simplify = F)
  g <- do.call(rbind, g)
  
  # terms <- c(terms, g$NAME[g$`FDR q-val` < 0.05])
  terms <- c(terms, g$NAME[(g$`FDR q-val` < 0.05) &
                             (abs(as.numeric(g$NES)) >= 2)])
}
rm(g); gc()


gsea_mat <- matrix(nrow = length(unique(terms)), ncol = length(ALL_TISSUES_GSEA),
                   dimnames = list(unique(terms), names(ALL_TISSUES_GSEA)))
gsea_pval_mat <- matrix(nrow = length(unique(terms)), ncol = length(ALL_TISSUES_GSEA),
                        dimnames = list(unique(terms), names(ALL_TISSUES_GSEA)))
for(t in unique(sapply(files, function(x) strsplit(strsplit(x,"_")[[1]][1],"\\.")[[1]][1]))){
  
  g <- sapply(files[grep(t,files)],
              function(f) read_delim(paste0(GSEA_DIR,"/",TYPE,"/",f),
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE), simplify = F)
  g <- do.call(rbind, g)
  
  gsea_mat[,
           names(ALL_TISSUES_GSEA)[ALL_TISSUES_GSEA == t]] <- as.numeric(g[match(rownames(gsea_mat),
                                                                                   g$NAME),]$NES)
  gsea_pval_mat[,
           names(ALL_TISSUES_GSEA)[ALL_TISSUES_GSEA == t]] <- as.numeric(g[match(rownames(gsea_pval_mat),
                                                                                   g$NAME),]$`FDR q-val`)
}


colnames(gsea_mat) <- tissue_labels_dash_all[colnames(gsea_mat)]

pdf(paste0(PLOT_DIR,"/GSEA_",
           TYPE,"_predictability_heatmap.pdf"), height = 7, width = 13)
ht <- Heatmap(gsea_mat, colorRamp2(breaks = c(-max(abs(gsea_mat), na.rm = T),
                                        0, max(abs(gsea_mat), na.rm = T)),
                             # colors = c("#034780","black","#bf0a26"),
                             # colors = c("#034780","white","#bf0a26"),
                             colors = c("#FF6347","white","#40547c"),
                             reverse = F),
        heatmap_legend_param = list(title = "Enrichment score"),
        # cluster_rows = FALSE,
        width = unit(10,"cm"), height = unit(10,"cm"),
        # row_title = "Modules",
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 12), row_dend_side = "right",
        row_names_side = "left", column_names_rot = 40)
draw(ht, padding = unit(c(0,10,0,0),"cm"))
dev.off()


stars.pval(gsea_pval_mat)


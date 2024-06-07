# Plot heatmap of slopes for all tissues

# SLOPE_DIR = "Outputs/5_Predictability/Age/Age_MaxSubset"
SLOPE_DIR = "Outputs/5_Predictability/Age"
# PLOT_DIR = "Plots/5_Predictability/Age/Age_MaxSubset"
PLOT_DIR = "Plots/5_Predictability/Age"
COR_THRE = "well_predicted"
TISSUES = "all"
INCL_70 = T

library(ComplexHeatmap)
library(circlize)
library(viridis)
source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")

if(INCL_70){
  slope_files <- list.files(SLOPE_DIR,
                            pattern = paste0("ageslope_",COR_THRE,".rds"),
                            full.names = T)
} else{
  slope_files <- list.files(SLOPE_DIR,
                            pattern = paste0("ageslope_",COR_THRE,"_no70.rds"),
                            full.names = T)
}
slopes <- sapply(slope_files, ReadRDS, simplify = F)
names(slopes) <- sapply(names(slopes),
                        function(x) strsplit(tail(strsplit(x, "/")[[1]],1),
                                             "_")[[1]][1])
rm(slope_files); gc()

# Restrict to tissues with signal
if(TISSUES != "all")
  slopes <- slopes[names(slopes) %in% MAIN_TISSUES]
common <- names(which(table(unlist(lapply(slopes,
                                          function(m) rownames(m)))) ==
                        length(slopes)))
slopes <- do.call(cbind, lapply(slopes, function(m) m[common,"Slope"]))

slopes <- slopes[order(rowMeans(slopes)),]

# Plot heatmap
# h <- Heatmap(slopes, col = colorRamp2(breaks = c(-0.012,0,0.012),
#                                       hcl_palette = "cividis", reverse = T),
#              row_title = "Genes", cluster_rows = F,
#              column_labels = tissue_labels_dash[colnames(slopes)],
#              show_row_names = FALSE, show_row_dend = FALSE,
#              heatmap_legend_param = list(title = "Slope",
#                                          direction = "horizontal",
#                                          title_position = "leftcenter",
#                                          legend_width = unit(4,"cm"),
#                                          title_gp = gpar(fontsize = 18),
#                                          labels_gp = gpar(fontsize = 16)),
#              # column_title = "Age-related predictability changes",
#              column_names_gp = gpar(fontsize = 16),
#              row_title_gp = gpar(fontsize = 20))
h <- Heatmap(t(slopes), col = colorRamp2(breaks = c(-0.012,0,0.012),
                                         # hcl_palette = "cividis", reverse = T),
                                         colors = c(age_palette["60-69"],"white",age_palette["20-29"])),
             column_title = "Genes", cluster_columns = F,
             row_labels = tissue_labels_dash[colnames(slopes)],
             show_column_names = FALSE, show_column_dend = FALSE,
             heatmap_legend_param = list(title = "Slope",
                                         direction = "horizontal",
                                         title_position = "leftcenter",
                                         legend_width = unit(4,"cm"),
                                         title_gp = gpar(fontsize = 18),
                                         labels_gp = gpar(fontsize = 16)),
             # column_title = "Age-related predictability changes",
             row_names_gp = gpar(fontsize = 16),
             column_title_gp = gpar(fontsize = 20))

if(INCL_70){
  pdf(paste0(PLOT_DIR,"/ageslope_comparison_maintissues.pdf"), height = 10)
} else{
  pdf(paste0(PLOT_DIR,"/ageslope_comparison_maintissues_no70.pdf"), height = 10)
}
h
dev.off()

if(INCL_70){
  if(TISSUES == "all"){
    WriteRDS(h, paste0(SLOPE_DIR,"/slope_heatmap_all.rds"))
  } else{
    WriteRDS(h, paste0(SLOPE_DIR,"/slope_heatmap_main.rds"))
  }
} else{
  WriteRDS(h, paste0("/slope_heatmap_main_no70.rds"))
}


# Heatmap of expression slopes
DE_files <- list.files("Outputs/3_GTExDataPrep/Differential_Expression",
                       pattern = "_DE")
DE <- sapply(paste0("Outputs/3_GTExDataPrep/Differential_Expression/",DE_files),
             ReadRDS, simplify = F)
names(DE) <- sapply(DE_files, function(x) strsplit(x,"_")[[1]][1])
DE <- lapply(DE, function(fit)
  DataENSGToSymbol(limma::topTable(fit, coef = "AGE",
                                   number = nrow(fit)),
                   remove_dup = T))
DE <- do.call(cbind, lapply(DE, function(m) m[common,"logFC"]))
rownames(DE) <- common

DE <- DE[rownames(slopes),colnames(slopes)[column_order(h)]]

pdf(paste0(PLOT_DIR,"/ageslope_comparison_DE_maintissues.pdf"), height = 10)
Heatmap(DE, col = colorRamp2(breaks = c(-0.013,0,0.013),
                             hcl_palette = "cividis", reverse = T),
        row_title = "Genes", cluster_rows = F, cluster_columns = F,
        show_row_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title = "Age slopes"),
        column_title = "Age-related expression changes",
        column_title_gp = gpar(fontsize = 20))
dev.off()


# Heatmap of variance slopes
var_files <- list.files("Outputs/3_GTExDataPrep/Differential_Expression",
                        pattern = "_var_regression")
var <- sapply(paste0("Outputs/3_GTExDataPrep/Differential_Expression/",var_files),
              ReadRDS, simplify = F)
names(var) <- sapply(var_files, function(x) strsplit(x,"_")[[1]][1])
var <- lapply(var, DataENSGToSymbol, remove_dup = T)
var <- do.call(cbind, lapply(var, function(m) m[common,"Slope"]))
rownames(var) <- common

var <- var[rownames(slopes),colnames(slopes)[column_order(h)]]

pdf(paste0(PLOT_DIR,"/ageslope_comparison_var_maintissues.pdf"), height = 10)
Heatmap(var, col = colorRamp2(breaks = c(-0.04,0,0.04),
                              hcl_palette = "cividis", reverse = T),
        row_title = "Genes", cluster_rows = F, cluster_columns = F,
        show_row_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title = "Age slopes"),
        column_title = "Age-related variance changes",
        column_title_gp = gpar(fontsize = 20))
dev.off()


# GO enrichment on common hits
GO <- list()
GO[["BP"]] <- GetGOEnrich(foreground = common,
                          background = na.omit(VectorENSGToSymbol(rownames(ReadRDS(
                            "Outputs/3_GTExDataPrep/Subset_Data/Adipose-Subcutaneous_20-29_sampled_data.rds")))),
                          go = "BP")
GO[["MF"]] <- GetGOEnrich(foreground = common,
                          background = na.omit(VectorENSGToSymbol(rownames(ReadRDS(
                            "Outputs/3_GTExDataPrep/Subset_Data/Adipose-Subcutaneous_20-29_sampled_data.rds")))),
                          go = "MF")
GO[["CC"]] <- GetGOEnrich(foreground = common,
                          background = na.omit(VectorENSGToSymbol(rownames(ReadRDS(
                            "Outputs/3_GTExDataPrep/Subset_Data/Adipose-Subcutaneous_20-29_sampled_data.rds")))),
                          go = "CC")
if(INCL_70){
  WriteRDS(GO,"Outputs/5_Predictability/Age/GOenrich_common_maintissues.rds")
} else{
  WriteRDS(GO,"Outputs/5_Predictability/Age/GOenrich_common_maintissues_no70.rds")
}



# Rowise correlations to correlation slopes
correlations <- list("Expression" = t(sapply(rownames(slopes),
                                             function(g) unlist(cor.test(x = slopes[g,],
                                                                         y = DE[g,],
                                                                         method = "pearson")[c("estimate",
                                                                                               "p.value")]))),
                     "Variance" = t(sapply(rownames(slopes),
                                           function(g) unlist(cor.test(x = slopes[g,],
                                                                       y = var[g,],
                                                                       method = "pearson")[c("estimate",
                                                                                             "p.value")]))))

expr_cor <- names(which(correlations$Expression[,"p.value"] < 0.05))
var_cor <- names(which(correlations$Variance[,"p.value"] < 0.05))
# slopes suck for all these genes



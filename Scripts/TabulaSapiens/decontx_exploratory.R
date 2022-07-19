# Ambient RNA decontamination - compare to uncorrected

.libPaths("../../Resources/Rlibs/R-4.0.3")
library(ggplot2)
library(ggpubr)
library(scran)
library(scuttle)

TISSUE = "Lung"

sce_10X <- readRDS(paste0("../../Outputs/Tabula_Sapiens/pertissue/10X/",TISSUE,".rds"))
sce_SS2 <- readRDS(paste0("../../Outputs/Tabula_Sapiens/pertissue/smartseq2/",TISSUE,".rds"))

cell_types <- intersect(names(which(table(sce_10X@colData$cell_ontology_class) > 100)),
                        names(which(table(sce_SS2@colData$cell_ontology_class) > 100)))

plot.data <- data.frame()
cor.data <- data.frame()

for(cell_type in cell_types){
  
  # restrict to cell type
  tmp_10X <- sce_10X[,sce_10X@colData$cell_ontology_class == cell_type]
  tmp_SS2 <- sce_SS2[,sce_SS2@colData$cell_ontology_class == cell_type]
  
  
  random_genes <- sample(1:nrow(tmp_10X), 2000)
  random_cells <- sample(1:ncol(tmp_10X), 100)
  
  plot.data <- rbind.data.frame(plot.data,
                                data.frame("raw" = as.numeric(assay(tmp_10X,
                                                                    "raw_counts")[random_genes,random_cells]),
                                           "decontx" = as.numeric(assay(tmp_10X,
                                                                        "decontXcounts")[random_genes,random_cells]),
                                           "Method" = "10X", "CellType" = cell_type))
  
  # library size correction
  tpm_raw <- sweep(as.matrix(assay(tmp_10X, "raw_counts")), 2,
                   STATS = colSums(as.matrix(assay(tmp_10X, "raw_counts")))/(10^6), FUN = "/")
  tpm_decontx <- sweep(as.matrix(assay(tmp_10X, "decontXcounts")), 2,
                       STATS = colSums(as.matrix(assay(tmp_10X, "decontXcounts")))/(10^6), FUN = "/")
  
  # library size correction (scran)
  sf_raw <- calculateSumFactors(tmp_10X, assay.type = "raw_counts")
  # sf_decontx <- calculateSumFactors(tmp_10X, assay.type = "decontXcounts")
  # log transformation and normalization
  tmp_10X <- logNormCounts(tmp_10X, size.factors = sf_raw, assay.type = "raw_counts",
                           name = "log_raw_counts")
  # tmp_10X <- logNormCounts(tmp_10X, size.factors = sf_decontx, assay.type = "decontXcounts",
  #                          name = "log_decontXcounts")
  rm(sf_raw,sf_decontx); gc()
  # identification of highly variable genes (in the raw counts)
  dec_raw <- modelGeneVar(tmp_10X, assay.type = "log_raw_counts")
  topvar_genes <- getTopHVGs(dec_raw, n = 5000)
  rm(dec_raw); gc()
  # identification of highly expressed genes (in the raw counts)
  means <- sort(rowMeans(tpm_raw), decreasing = T)
  top_genes <- head(intersect(head(names(means),5000), topvar_genes),2000)
  rm(means); gc()
  
  # cormat_raw <- cor(t(as.matrix(assay(tmp_10X, "log_raw_counts"))[topvar_genes,random_cells]))
  cormat_raw <- cor(t(tpm_raw[top_genes,random_cells]))
  # cormat_decontx <- cor(t(as.matrix(assay(tmp_10X, "log_decontXcounts"))[topvar_genes,random_cells]))
  cormat_decontx <- cor(t(tpm_decontx[top_genes,random_cells]))
  cor.data <- rbind.data.frame(cor.data,
                               data.frame("raw" = cormat_raw[upper.tri(cormat_raw)],
                                          "decontx" = cormat_decontx[upper.tri(cormat_decontx)],
                                          "Method" = "10X", "CellType" = cell_type))
  
  rm(cormat_raw,cormat_decontx,random_cells,topvar_genes); gc()
  
  
  random_cells <- sample(1:ncol(tmp_SS2), 100)
  
  plot.data <- rbind.data.frame(plot.data,
                                data.frame("raw" = as.numeric(assay(tmp_SS2,
                                                                    "raw_counts")[random_genes,random_cells]),
                                           "decontx" = as.numeric(assay(tmp_SS2,
                                                                        "decontXcounts")[random_genes,random_cells]),
                                           "Method" = "Smart-seq2", "CellType" = cell_type))
  
  # library size correction
  tpm_raw <- sweep(as.matrix(assay(tmp_SS2, "raw_counts")), 2,
                   STATS = colSums(as.matrix(assay(tmp_SS2, "raw_counts")))/(10^6), FUN = "/")
  tpm_decontx <- sweep(as.matrix(assay(tmp_SS2, "decontXcounts")), 2,
                       STATS = colSums(as.matrix(assay(tmp_SS2, "decontXcounts")))/(10^6), FUN = "/")
  # library size correction (scran)
  sf_raw <- calculateSumFactors(tmp_SS2, assay.type = "raw_counts")
  # sf_decontx <- calculateSumFactors(tmp_SS2, assay.type = "decontXcounts")
  # log transformation and normalization
  tmp_SS2 <- logNormCounts(tmp_SS2, size.factors = sf_raw, assay.type = "raw_counts",
                           name = "log_raw_counts")
  # tmp_SS2 <- logNormCounts(tmp_SS2, size.factors = sf_decontx, assay.type = "decontXcounts",
  #                          name = "log_decontXcounts")
  rm(sf_raw,sf_decontx); gc()
  # identification of highly variable genes (in the raw counts)
  dec_raw <- modelGeneVar(tmp_SS2, assay.type = "log_raw_counts")
  topvar_genes <- getTopHVGs(dec_raw, n = 5000)
  rm(dec_raw); gc()
  # identification of highly expressed genes (in the raw counts)
  means <- sort(rowMeans(tpm_raw), decreasing = T)
  top_genes <- head(intersect(head(names(means),5000), topvar_genes),2000)
  rm(means); gc()
  
  # cormat_raw <- cor(t(as.matrix(assay(tmp_SS2, "log_raw_counts"))[topvar_genes,random_cells]))
  cormat_raw <- cor(t(tpm_raw[top_genes,random_cells]))
  # cormat_decontx <- cor(t(as.matrix(assay(tmp_SS2, "log_decontXcounts"))[topvar_genes,random_cells]))
  cormat_decontx <- cor(t(tpm_decontx[top_genes,random_cells]))
  cor.data <- rbind.data.frame(cor.data,
                               data.frame("raw" = cormat_raw[upper.tri(cormat_raw)],
                                          "decontx" = cormat_decontx[upper.tri(cormat_decontx)],
                                          "Method" = "Smart-seq2", "CellType" = cell_type))
  
  rm(cormat_decontx,cormat_raw,random_cells,topvar_genes); gc()
  rm(tpm_raw,tpm_decontx); gc()
  rm(tmp_10X,tmp_SS2); gc()
}


# plot scatterplots corrected vs uncorrected

plots <- list()
plots[[1]] <- ggplot(subset(plot.data, Method == "10X")) +
  geom_point(aes(x = log2(1+raw), y = log2(1+decontx))) +
  xlab("Raw counts - log2-transformed") +
  ylab("Corrected counts - log2-transformed") +
  facet_wrap(~ CellType, scales = "free") +
  ggtitle(TISSUE," - 10X") + theme(text = element_text(size = 20))
plots[[2]] <- ggplot(subset(plot.data, Method == "Smart-seq2")) +
  geom_point(aes(x = log2(1+raw), y = log2(1+decontx))) +
  xlab("Raw counts - log2-transformed") +
  ylab("Corrected counts - log2-transformed") +
  facet_wrap(~ CellType, scales = "free") +
  ggtitle(TISSUE," - Smart-seq2") + theme(text = element_text(size = 20))


png(paste0("../../Plots/Tabula_Sapiens/QC/",TISSUE,"_decontx.png"),
    # width = 900,
    width = 500,
    height = 800)
ggarrange(plotlist = plots, nrow = 2)
dev.off()

rm(plots, plot.data); gc()


# plot gene-gene correlations

plots <- list()
plots[[1]] <- ggplot(subset(cor.data, Method == "10X")) +
  geom_density(aes(x = raw, color = "Raw counts"), size = 1) +
  geom_density(aes(x = decontx, color = "Corrected counts"), size = 1) +
  xlab("Gene-gene correlations") + coord_cartesian(xlim = c(-1,1)) +
  facet_wrap(~ CellType, scales = "free_y") +
  ggtitle(TISSUE," - 10X") + theme(text = element_text(size = 20),
                                   legend.title = element_blank())
plots[[2]] <- ggplot(subset(cor.data, Method == "Smart-seq2")) +
  geom_density(aes(x = raw, color = "Raw counts"), size = 1) +
  geom_density(aes(x = decontx, color = "Corrected counts"), size = 1) +
  xlab("Gene-gene correlations") + coord_cartesian(xlim = c(-1,1)) +
  facet_wrap(~ CellType, scales = "free_y") + 
  ggtitle(TISSUE," - Smart-seq2") + theme(text = element_text(size = 20),
                                          legend.title = element_blank())

png(paste0("../../Plots/Tabula_Sapiens/QC/",TISSUE,"_decontx_genecor_topvargenes.png"),
    width = 700,
    # width = 950,
    # height = 600,
    height = 800)
ggarrange(plotlist = plots, nrow = 2, common.legend = T, legend = "right")
dev.off()

rm(cor.data,plots,cell_type,cell_types); gc()


# cluster cells before and after decontX (PCA)
library(scater)
plots <- list()

# raw counts
sce_SS2 <- logNormCounts(sce_SS2, assay.type = "raw_counts")
sce_SS2 <- runPCA(sce_SS2, exprs_values = "logcounts", name = "PCA")
plots[[1]] <- plotPCA(sce_SS2, colour_by = "cell_ontology_class",
                      shape_by = "donor") + ggtitle("Smart-seq2 raw counts") +
  guides(shape = F)
# corrected counts
sce_SS2 <- logNormCounts(sce_SS2, assay.type = "decontXcounts",
                         name = "decontXlogcounts")
sce_SS2 <- runPCA(sce_SS2, exprs_values = "decontXlogcounts", name = "PCA")
plots[[2]] <- plotPCA(sce_SS2, colour_by = "cell_ontology_class",
                      shape_by = "donor") + ggtitle("Smart-seq2 decontX") +
  guides(shape = F)

# raw counts
sce_10X <- logNormCounts(sce_10X, assay.type = "raw_counts")
sce_10X <- runPCA(sce_10X, exprs_values = "logcounts", name = "PCA")
plots[[3]] <- plotPCA(sce_10X, colour_by = "cell_ontology_class",
                      shape_by = "donor") + ggtitle("10X raw counts") +
  guides(shape = F)
# corrected counts
sce_10X <- logNormCounts(sce_10X, assay.type = "decontXcounts",
                         name = "decontXlogcounts")
sce_10X <- runPCA(sce_10X, exprs_values = "decontXlogcounts", name = "PCA")
plots[[4]] <- plotPCA(sce_10X, colour_by = "cell_ontology_class",
                      shape_by = "donor") + ggtitle("10X decontX") +
  guides(shape = F)

pdf(paste0("../../Plots/Tabula_Sapiens/QC/",TISSUE,"_decontx_pca.pdf"),
    width = 10, height = 9)
ggarrange(ggarrange(plotlist = plots[1:2], nrow = 1, ncol = 2, align = "h", common.legend = T, legend = "bottom"),
          ggarrange(plotlist = plots[3:4], nrow = 1, ncol = 2, align = "h", common.legend = T, legend = "bottom"),
          nrow = 2, align = "v")
dev.off()


rm(plots, sce_10X, sce_SS2); gc()

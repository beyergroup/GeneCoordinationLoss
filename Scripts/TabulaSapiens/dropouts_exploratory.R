# If with Smart-seq we have a higher saturation (sequence more deeply but
# detect the same number of genes), we expect dropouts to be biological and
# thus more stable across individual cells.

library(ggplot2)

TISSUE = "Vasculature"

sce_10X <- readRDS(paste0("../../Outputs/Tabula_Sapiens/pertissue/10X/",TISSUE,".rds"))
sce_SS2 <- readRDS(paste0("../../Outputs/Tabula_Sapiens/pertissue/smartseq2/",TISSUE,".rds"))

# # remove genes that dropout across all cells of tissue
# sce_10X <- sce_10X[rowSums(assay(sce_10X,"raw_counts") != 0) > 0,]
# sce_SS2 <- sce_SS2[rowSums(assay(sce_SS2,"raw_counts") != 0) > 0,]

# # remove genes that never dropout in any cell of tissue
# sce_10X <- sce_10X[rowSums(assay(sce_10X,"raw_counts") == 0) > 0,]
# sce_SS2 <- sce_SS2[rowSums(assay(sce_SS2,"raw_counts") == 0) > 0,]

# look at abundant cell types in both datasets
cell_types <- intersect(names(which(table(sce_10X@colData$cell_ontology_class) > 100)),
                        names(which(table(sce_SS2@colData$cell_ontology_class) > 100)))

plot.data <- data.frame()
drop.count <- data.frame()

for(cell_type in cell_types){
  
  # restrict to cell type
  tmp_10X <- sce_10X[,sce_10X@colData$cell_ontology_class == cell_type]
  tmp_SS2 <- sce_SS2[,sce_SS2@colData$cell_ontology_class == cell_type]
  
  # subsample smallest dataset to make distributions comparable
  ncells <- min(c(ncol(tmp_10X),ncol(tmp_SS2)))
  tmp_10X <- tmp_10X[,sample(1:ncol(tmp_10X),ncells)]
  tmp_SS2 <- tmp_SS2[,sample(1:ncol(tmp_SS2),ncells)]
  
  drop.count <- rbind.data.frame(drop.count,
                                 data.frame("GlobalDropouts" = sum(assay(tmp_10X, "raw_counts") == 0)/(ncol(tmp_10X)*nrow(tmp_10X)),
                                            "Method" = "10X", "CellType" = cell_type, "CellCount" = ncol(tmp_10X)))
  drop.count <- rbind.data.frame(drop.count,
                                 data.frame("GlobalDropouts" = sum(assay(tmp_SS2, "raw_counts") == 0)/(ncol(tmp_SS2)*nrow(tmp_SS2)),
                                            "Method" = "Smart-seq2", "CellType" = cell_type,
                                            "CellCount" = ncol(tmp_SS2)))
  
  # remove genes that dropout across all cells of tissue and restrict to cell type
  tmp_10X <- tmp_10X[rowSums(assay(sce_10X,"raw_counts") != 0) > 0,]
  tmp_SS2 <- tmp_SS2[rowSums(assay(sce_SS2,"raw_counts") != 0) > 0,]
  
  
  plot.data <- rbind.data.frame(plot.data,
                                data.frame("Fraction" = rowSums(assay(tmp_10X,"raw_counts") == 0)/ncol(tmp_10X),
                                           "Method" = "10X", "CellType" = cell_type))
  plot.data <- rbind.data.frame(plot.data,
                                data.frame("Fraction" = rowSums(assay(tmp_SS2,"raw_counts") == 0)/ncol(tmp_SS2),
                                           "Method" = "Smart-seq2", "CellType" = cell_type))
  
  rm(tmp_10X,tmp_SS2); gc()
}

pdf(paste0("../../Plots/Tabula_Sapiens/QC/",TISSUE,"_dropouts_subsampled.pdf"), width = 9)

ggplot(plot.data) +
  geom_violin(aes(y = Fraction, x = Method, fill = Method)) +
  coord_cartesian(ylim = c(0.75,1)) +
  ylab("Fraction of dropout cells") +
  facet_wrap(~ CellType, scales = "free_y") +
  guides(fill = FALSE) +
  ggtitle(paste0(TISSUE," - genes detected in tissue")) +
  theme(text = element_text(size = 20))

ggplot(drop.count) +
  geom_bar(aes(x = CellType, y = GlobalDropouts, fill = Method),
           stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(x = CellType, y = (GlobalDropouts+0.01), label = paste0("N = ",CellCount), group = Method),
            position = position_dodge(width = .9), check_overlap = TRUE, hjust = 0) +
  ylab("Data sparsity") + coord_flip() +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
                     limits = c(0,1.15)) +
  ggtitle(paste0(TISSUE," - whole data")) +
  theme(text = element_text(size = 20),
        legend.position = "bottom")
dev.off()

rm(plot.data,drop.count,sce_10X,sce_SS2); gc()

# Normalization strategies

.libPaths("../../Resources/Rlibs/R-4.0.3")
library(SingleCellExperiment)
library(scuttle)
library(ggplot2)
library(ggpubr)
library(reshape2)


TISSUE = "Pancreas"

sce_10X <- readRDS(paste0("../../Outputs/Tabula_Sapiens/CellFiltered/10X_",TISSUE,".rds"))
sce_SS2 <- readRDS(paste0("../../Outputs/Tabula_Sapiens/CellFiltered/smartseq2_",TISSUE,".rds"))
assay(sce_10X,"raw_counts") <- NULL
assay(sce_SS2,"raw_counts") <- NULL
assay(sce_10X,"X") <- NULL
assay(sce_SS2,"X") <- NULL
gc()


# evaluate normalization within cell type
cell_types <- intersect(names(which(table(sce_10X@colData$cell_ontology_class) > 100)),
                        names(which(table(sce_SS2@colData$cell_ontology_class) > 100)))

plot.data <- data.frame()
cor.data <- data.frame()

# pdf(paste0("../../Plots/Tabula_Sapiens/QC/",TISSUE,"_raw_genequant.pdf"),
#     height = 12, width = 10)

for(cell_type in cell_types){
  
  # plot.lines <- data.frame()
  
  # restrict to cell type
  tmp_10X <- sce_10X[,sce_10X@colData$cell_ontology_class == cell_type]
  tmp_SS2 <- sce_SS2[,sce_SS2@colData$cell_ontology_class == cell_type]
  
  # expression profiles of random cells
  cells <- c(sample(head(order(tmp_10X@colData$n_counts_UMIs, decreasing = T),10),1), # high libsize
             sample(head(order(tmp_10X@colData$n_counts_UMIs, decreasing = F),10),1)) # low libsize
  plot.data <- rbind.data.frame(plot.data,
                                data.frame("Cell1" = as.matrix(assay(tmp_10X,"decontXcounts"))[,cells[1]],
                                           "Cell2" = as.matrix(assay(tmp_10X,"decontXcounts"))[,cells[2]],
                                           "CellType" = cell_type,
                                           "Method" = "10X"))
  cells <- c(sample(head(order(tmp_SS2@colData$n_counts_UMIs, decreasing = T),10),1), # high libsize
             sample(head(order(tmp_SS2@colData$n_counts_UMIs, decreasing = F),10),1)) # low libsize
  plot.data <- rbind.data.frame(plot.data,
                                data.frame("Cell1" = as.matrix(assay(tmp_SS2,"decontXcounts"))[,cells[1]],
                                           "Cell2" = as.matrix(assay(tmp_SS2,"decontXcounts"))[,cells[2]],
                                           "CellType" = cell_type,
                                           "Method" = "Smart-seq2"))
  
  # expression quantiles (5 quantiles)
  means <- rowMeans(log2(1+as.matrix(assay(tmp_10X,"decontXcounts"))))
  tmp_10X <- tmp_10X[means != 0,]; means <- means[means != 0]
  quantiles <- c(0,quantile(means, prob = c(0.2,0.4,0.6,0.8)),max(means))
  q <- sapply(means, function(x) sum(x > quantiles))
  rm(means,quantiles); gc()
  
  for(qq in 1:5){
    genes <- names(which(q == qq))
    cor <- apply(log2(1+as.matrix(assay(tmp_10X,"decontXcounts")))[genes,], 1,
                 function(xx) cor(x = log10(tmp_10X@colData$n_counts_UMIs),
                                  y = xx,
                                  method = "spearman"))
    cor.data <- rbind.data.frame(cor.data,
                                 data.frame("Cor" = cor, "GeneGroup" = qq,
                                            "CellType" = cell_type, "Method" = "10X"))
    rm(genes,cor); gc()
    # # fit kernel smooth function per gene (expression ~ log10 lib size)
    # fit <- apply(log2(1+as.matrix(assay(tmp_10X,"decontXcounts")))[genes,], 1,
    #              function(xx) ksmooth(x = log10(tmp_10X@colData$n_counts_UMIs), y = xx,
    #                                  bandwidth = 20*bw.SJ(log10(tmp_10X@colData$n_counts_UMIs))))
    # # median and IQR of lines (across genes in group)
    # plot.lines <- rbind.data.frame(plot.lines,
    #                                data.frame("Mean" = rowMedians(do.call(cbind, lapply(fit, function(l) l$y))),
    #                                           "IQR1" = rowQuantiles(do.call(cbind,
    #                                                                         lapply(fit, function(l) l$y)),
    #                                                                 probs = 0.25),
    #                                           "IQR3" = rowQuantiles(do.call(cbind,
    #                                                                         lapply(fit, function(l) l$y)),
    #                                                                 probs = 0.75),
    #                                           "LibSize" = rowMeans(do.call(cbind, lapply(fit, function(l) l$x))),
    #                                           "GeneGroup" = qq, "CellType" = cell_type,
    #                                           "Method" = "10X"))
    # rm(genes,fit); gc()
  }
  
  rm(q,tmp_10X); gc()
  
  
  means <- rowMeans(log2(1+as.matrix(assay(tmp_SS2,"decontXcounts"))))
  tmp_SS2 <- tmp_SS2[means != 0,]; means <- means[means != 0]
  quantiles <- c(0,quantile(means, prob = c(0.2,0.4,0.6,0.8)),max(means))
  q <- sapply(means, function(x) sum(x > quantiles))
  rm(means,quantiles); gc()
  
  for(qq in 1:5){
    genes <- names(which(q == qq))
    cor <- apply(log2(1+as.matrix(assay(tmp_SS2,"decontXcounts")))[genes,], 1,
                 function(xx) cor(x = log10(tmp_SS2@colData$n_counts_UMIs),
                                  y = xx,
                                  method = "spearman"))
    cor.data <- rbind.data.frame(cor.data,
                                 data.frame("Cor" = cor, "GeneGroup" = qq,
                                            "CellType" = cell_type, "Method" = "Smart-seq2"))
    rm(genes,cor); gc()
    # # fit kernel smooth function per gene (expression ~ log10 lib size)
    # fit <- apply(log2(1+as.matrix(assay(tmp_SS2,"decontXcounts")))[genes,], 1,
    #              function(xx) ksmooth(x = log10(tmp_SS2@colData$n_counts_UMIs), y = xx,
    #                                   bandwidth = 20*bw.SJ(log10(tmp_SS2@colData$n_counts_UMIs))))
    # # median and IQR of lines (across genes in group)
    # plot.lines <- rbind.data.frame(plot.lines,
    #                                data.frame("Mean" = rowMedians(do.call(cbind, lapply(fit, function(l) l$y))),
    #                                           "IQR1" = rowQuantiles(do.call(cbind,
    #                                                                         lapply(fit, function(l) l$y)),
    #                                                                 probs = 0.25),
    #                                           "IQR3" = rowQuantiles(do.call(cbind,
    #                                                                         lapply(fit, function(l) l$y)),
    #                                                                 probs = 0.75),
    #                                           "LibSize" = rowMeans(do.call(cbind, lapply(fit, function(l) l$x))),
    #                                           "GeneGroup" = as.character(qq), "CellType" = cell_type,
    #                                           "Method" = "Smart-seq2"))
    # rm(genes,fit); gc()
  }
  
  rm(q,tmp_SS2); gc()
  
  # # Generate plot of trend expression ~ lib size
  # plots <- list()
  # plots[[1]] <- ggplot(subset(plot.lines, Method == "10X")) +
  #   facet_wrap(~ GeneGroup, scales = "free_y", ncol = 1,
  #              labeller = as_labeller(function(x) paste0("Q",x))) +
  #   geom_ribbon(aes(ymin  = IQR1, ymax = IQR3, x = LibSize, fill = GeneGroup)) +
  #   geom_line(aes(y = Mean, x = LibSize, group = GeneGroup), size = 2) +
  #   ggtitle(paste(TISSUE,cell_type)," - 10X") +
  #   xlab("Library size (log scale)") +
  #   ylab("Log counts") + guides(fill = FALSE) +
  #   scale_fill_viridis_d(begin = .2, end = .9) +
  #   theme(text = element_text(size = 20))
  # plots[[2]] <- ggplot(subset(plot.lines, Method == "Smart-seq2")) +
  #   facet_wrap(~ GeneGroup, scales = "free_y", ncol = 1,
  #              labeller = as_labeller(function(x) paste0("Q",x))) +
  #   geom_ribbon(aes(ymin  = IQR1, ymax = IQR3, x = LibSize, fill = GeneGroup)) +
  #   geom_line(aes(y = Mean, x = LibSize, group = GeneGroup), size = 2) +
  #   ggtitle("","- Smart-seq2") +
  #   xlab("Library size (log scale)") +
  #   ylab("Log counts") + guides(fill = FALSE) +
  #   scale_fill_viridis_d(begin = .2, end = .9) +
  #   theme(text = element_text(size = 20),
  #         plot.title = element_blank(),
  #         plot.background = element_blank())
  # 
  # print(ggarrange(plotlist = plots, ncol = 2, common.legend = T,
  #                 legend = "right", align = "h"))
  
}

# dev.off()

plot.data <- subset(plot.data, (Cell1 != 0) | (Cell2 != 0))

# Generate plot of cell vs cell
plots <- list()
plots[[1]] <- ggplot(subset(plot.data, Method == "10X")) +
  facet_wrap(~ CellType) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(aes(x = Cell1, y = Cell2)) +
  xlab("High LibSize cell") +
  ylab("Low LibSize cell") +
  ggtitle(TISSUE," - 10X") +
  theme(text = element_text(size = 20))
plots[[2]] <- ggplot(subset(plot.data, Method == "Smart-seq2")) +
  facet_wrap(~ CellType) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(aes(x = Cell1, y = Cell2)) +
  xlab("High LibSize cell") +
  ylab("Low LibSize cell") +
  ggtitle(TISSUE," - Smart-seq2") +
  theme(text = element_text(size = 20))

png(paste0("../../Plots/Tabula_Sapiens/QC/",TISSUE,"_raw_cellvscell.png"),
    # width = 900, height = 700) # Vasculature
    # width = 600, height = 700) # Pancreas
    width = 600, height = 1000) # Lung
ggarrange(plotlist = plots, nrow = 2)
dev.off()

rm(plots); gc()


# Generate plot grid of correlation densities
plots <- list()
plots[[1]] <- ggplot(subset(cor.data, Method == "10X")) +
  geom_density(aes(x = Cor, fill = as.character(GeneGroup))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_viridis_d(begin = .2, end = .9) +
  xlab("Correlation coefficient") +
  ggtitle(TISSUE," - 10X") + guides(fill = F) +
  facet_grid(GeneGroup ~ CellType,
             labeller = labeller(GeneGroup = function(x) paste0("Q",x)),
             scales = "free") +
  coord_cartesian(xlim = c(-1,1)) +
  theme(text = element_text(size = 20))
plots[[2]] <- ggplot(subset(cor.data, Method == "Smart-seq2")) +
  geom_density(aes(x = Cor, fill = as.character(GeneGroup))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_viridis_d(begin = .2, end = .9) +
  xlab("Correlation coefficient") + guides(fill = F) +
  ggtitle(TISSUE," - Smart-seq2") +
  facet_grid(GeneGroup ~ CellType,
             labeller = labeller(GeneGroup = function(x) paste0("Q",x)),
             scales = "free") +
  coord_cartesian(xlim = c(-1,1)) +
  theme(text = element_text(size = 20))

pdf(paste0("../../Plots/Tabula_Sapiens/QC/",TISSUE,"_raw_corr.pdf"),
    # width = 13, # Lung
    width = 6.5, # Pancreas
    height = 15)
ggarrange(plotlist = plots, nrow = 2)
dev.off()

rm(plots); gc()
rm(plot.data,plot.lines,cor.data,sce_10X,sce_SS2); gc()

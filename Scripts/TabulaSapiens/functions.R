NormEval_ExprScatter <- function(sce, assay, title){
  
  library(ggplot2)
  
  cells <- c(sample(head(order(sce@colData$n_counts_UMIs, decreasing = T),10),1), # high libsize
             sample(head(order(sce@colData$n_counts_UMIs, decreasing = F),10),1)) # low libsize
  
  plot.data <- data.frame("Cell1" = as.matrix(assay(sce,assay))[,cells[1]],
                          "Cell2" = as.matrix(assay(sce,assay))[,cells[2]])
  
  plot <- ggplot(plot.data) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_point(aes(x = Cell1, y = Cell2)) +
    xlab("High LibSize cell") +
    ylab("Low LibSize cell") +
    ggtitle(title[1],title[2]) +
    theme(text = element_text(size = 20))
  
  return(plot)
}


NormEval_LibTrend <- function(sce, assay, raw_assay, Q, title){
  
  library(ggplot2)
  
  means <- rowMeans(log2(1+as.matrix(assay(sce,raw_assay))))
  quantiles <- c(0,quantile(means, prob = seq(from = 0, to = 1, length.out = (Q+1)))[-1])
  q <- sapply(means, function(x) sum(x > quantiles))
  rm(means,quantiles); gc()
  
  cor.data <- data.frame()
  for(qq in 1:Q){
    genes <- names(which(q == qq))
    cor <- apply(log2(1+as.matrix(assay(sce,assay)))[genes,], 1,
                 function(xx) cor(x = log10(sce@colData$n_counts_UMIs),
                                  y = xx,
                                  method = "spearman"))
    cor.data <- rbind.data.frame(cor.data,
                                 data.frame("Cor" = cor, "GeneGroup" = qq))
    rm(genes,cor); gc()
  }
  
  plot <- ggplot(cor.data) +
    geom_density(aes(x = Cor, fill = as.character(GeneGroup))) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_viridis_d(begin = .2, end = .9) +
    xlab("Correlation coefficient") +
    ggtitle(title[1],title[2]) + guides(fill = F) +
    facet_grid(GeneGroup ~ .,
               labeller = labeller(GeneGroup = function(x) paste0("Q",x)),
               scales = "free_y") +
    coord_cartesian(xlim = c(-1,1)) +
    theme(text = element_text(size = 20),
          axis.title.y = element_blank(),
          plot.title = element_text(color = "transparent"),
          plot.background = element_blank(),
          plot.subtitle = element_text(color = "transparent"))
  
  return(plot)
}


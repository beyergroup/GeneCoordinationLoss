# topGO analysis --------------------------------------------------------------

# foreground, background: character vectors with genes in foreground and background
# go: either "BP" (biological processes), "MF" (molecular functions), "CC" (cellular components)
GetGOEnrich <- function(foreground, background, go, algorithm = "weight01",
                        pval_cutoff = 0.01, enrich_cutoff = 0.5){
  
  library(topGO)
  genes <- as.numeric(background %in% foreground)
  names(genes) <- background
  
  sampleGOdata <- new(Class = "topGOdata", ontology = go, allGenes = genes,
                      geneSelectionFun = function(x){x == 1}, nodeSize = 20, annotationFun = annFUN.org,
                      mapping = "org.Hs.eg.db", ID = "symbol")
  resultFisher <- runTest(sampleGOdata, algorithm = algorithm, statistic = "fisher")
  allRes <- GenTable(sampleGOdata, pval = resultFisher, numChar = 250,
                     orderBy = "pval", ranksOf = "pval",
                     topNodes = length(usedGO(sampleGOdata)))
  allRes$log2Enrichment <- log2(allRes$Significant/allRes$Expected)
  
  allRes <- allRes[allRes$pval < pval_cutoff & abs(allRes$log2Enrichment) > enrich_cutoff,]
  allRes <- allRes[order(allRes$log2Enrichment, decreasing = F),]
  allRes$Term <- factor(as.character(allRes$Term), levels = as.character(allRes$Term))
  
  return(allRes)
}


# GOenrich: output from GetGOEnrich()
# col: color for bars
PlotGOEnrich <- function(GOenrich, col, title){
  
  library(ggplot2)
  
  p <- ggplot(GOenrich) +
    geom_bar(aes(x = Term, y = log2Enrichment, alpha = -log10(as.numeric(pval))),
             stat = "identity", fill = col) + guides(alpha = F) +
    ggtitle(title) + coord_flip() + theme_classic() +
    theme(text = element_text(size = 20))
  
  return(p)
}


# Network functions -----------------------------------------------------------

# Visualize groups of genes in the network
# net: igraph object, undirected
NetVis <- function(net, genes = NULL, col){
  
  E(net)$weight <- abs(E(net)$weight)
  
  # adjust visual details
  V(net)$size <- 2
  V(net)$frame.color <- "white"
  V(net)$color <- "orange"
  V(net)$label <- "" 
  E(net)$arrow.mode <- 0
  
  # highlight provided genes
  
  
  # pick layout
  # ws  <-  c(1, rep(100, ecount(net)-1))
  # lw <- layout_with_fr(net, weights=ws)
  
  plot(net, layout = layout_kk)
  
  return()
}
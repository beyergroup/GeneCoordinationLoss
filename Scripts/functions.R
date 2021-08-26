# Data handling ---------------------------------------------------------------

DataENSGToSymbol <- function(data,
  conversion_file = "Resources/ensembl_idversion_GTExDESeq2_symbolChrStart.txt",
  remove_dup = F){
  
  table <- read.delim(conversion_file)
  rownames(data) <- sapply(rownames(data),
                           function(c) strsplit(c, split = "\\.")[[1]][1])
  rnames <- table[match(rownames(data), table$ensembl_gene_id),"symbol"]
  data <- data[!is.na(rnames),]
  rnames <- rnames[!is.na(rnames)]
  
  if(remove_dup){
    data <- data[-which(rnames %in% rnames[duplicated(rnames)]),]
    rnames <- rnames[-which(rnames %in% rnames[duplicated(rnames)])]
  }
  
  rownames(data) <- rnames
  rm(table,rnames); gc()
  
  return(data)
}


VectorENSGToSymbol <- function(x,
  conversion_file = "Resources/ensembl_idversion_GTExDESeq2_symbolChrStart.txt"){
  
  table <- read.delim(conversion_file)
  x <- sapply(x, function(c) strsplit(c, split = "\\.")[[1]][1])
  symbols <- table[match(x, table$ensembl_gene_id),"symbol"]
  rm(table); gc()
  
  return(symbols)
}


# GO analysis -----------------------------------------------------------------

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


# Given a GO ID, returns the genes annotated to that term
GetGOGenes <- function(GOID){
  
  library(biomaRt)
  
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  genes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","go_id"),
                 filters = "go", values = GOID, mart = ensembl)
  
  return(genes)
}



# Network functions -----------------------------------------------------------


# Remove residual edges in edge pairs
# threshold: ratio of weights (min edge)/(edge pair)

RemoveResEdges <- function(net, threshold = 0.1){
  
  if(any(colnames(net) != rownames(net)))
    stop("Network rows and columns must match.")
  
  # indices of upper triangle of net adj matrix, to iterate over
  # each row corresponds to one edge [i,j] in an edge pair ([i,j], [j,i])
  ind <- which(upper.tri(net), arr.ind = T)
  
  # compute ratio of weights for each edge pair
  ratios <- apply(ind, 1,
      function(x) min(abs(c(net[x[1],x[2]],
                            net[x[2],x[1]])))/sum(abs(c(net[x[1],x[2]],
                                                        net[x[2],x[1]]))))
  
  # select edges with ratios below the threshold - residual edges
  to_remove <- ind[(ratios < threshold) & (!is.na(ratios)),]
  rm(ind, ratios); gc()
  to_remove <- t(apply(to_remove, 1, function(x)
    switch(which.min(c(abs(net[x[1],x[2]]),abs(net[x[2],x[1]]))),
           c(x[1],x[2]),
           c(x[2],x[1]))))
  
  # set residual edges to 0
  net[to_remove] <- 0
  
  return(net)
}


# Predict expression in data given network
PredictNet <- function(net, centered_data, maxiter){
  
  # restrict to network targets and predictors
  targets <- intersect(rownames(centered_data),rownames(net))
  predictors <- intersect(rownames(centered_data),colnames(net))
  net <- net[targets, predictors]
  
  # multiply by network coefficients
  prediction <- centered_data
  for(i in 1:maxiter){
    message("Iteration ",i)
    prediction[targets,] <- net[targets,predictors] %*% prediction[predictors,]
    i <- i+1
  }
  
  return(prediction[targets,])
}




PlotNet <- function(graph, title, layout = "fr"){
  
  p <- plot(current_graph, layout = switch(layout,
                                           "lgl" = layout_with_lgl,
                                           "fr" = layout_with_fr,
                                           "nicely" = layout_nicely),
       vertex.size = 5, edge.curved = .1, vertex.label.dist = 1,
       edge.arrow.size = .2, vertex.frame.color = NA,
       vertex.label.color = "grey50", vertex.label.family = "Helvetica",
       vertex.label.cex = .8, main = title, margin = 0)
  
  return(p)
}

age_palette <- c("20-29" = "#126782", "30-39" = "#219EBC",
                 "50-59" = "#FD9E02", "60-69" = "#FB8500")

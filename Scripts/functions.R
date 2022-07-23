# Clustering functions --------------------------------------------------------

euc_dist <- function(m) {
  mtm <- Matrix::tcrossprod(m)
  sq <- rowSums(m*m)
  return(sqrt(outer(sq,sq,"+") - 2*mtm))
}


# Predictability functions ----------------------------------------------------

ComputeMSE <- function(pred, obs){
  
  se <- (pred-obs)^2
  e <- mean(na.omit(se))
  
  return(e)
}


ComputeNMSE <- function(pred, obs){
  
  e <- ComputeMSE(pred, obs)
  ne <- e/(mean(obs)^2)
  
  return(ne)
}


ComputeNMSEPerSample <- function(pred, obs){
  
  se <- (pred-obs)^2
  nse <- se/(obs^2)
  
  return(nse)
}


ComputeRMSEPerSample <- function(pred, obs){
  
  se <- (pred-obs)^2
  rse <- sqrt(se)
  
  return(rse)
}

# condition: vector of "interesting" prediction errors (higher errors here
# result in positive LFCs)
# control: vector of control prediction errors

ComputeErrorLFC <- function(condition, control){
  
  return(2*(log(condition) - log(control))/(log(condition) + log(control)))
  
}


# condition: vector of "interesting" prediction correlations (higher cors here
# result in negative LFCs)
# control: vector of control prediction correlations

ComputeCorrelationLFC <- function(condition, control){
  
  condition <- 1-condition
  control <- 1-control
  return(2*(log(condition) - log(control))/(log(condition) + log(control)))
}


ComputeCorrelationLFC_adj <- function(condition, control){
  
  condition <- 2-condition
  control <- 2-control
  return((log2(condition) - log2(control)) / 
           ((log2(condition) + log2(control)) / 2))
}



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


# Based on parameters wanted for network, return file according to personal
# naming conventions

GenerateNetworkName <- function(param){
  
  name <- paste0("Adjacency_",
                 c("undirected","directed")[param[["Directed"]]+1],
                 "_weights",param[["Weights"]],
                 c("","_rownormalized")[param[["Normalized"]]+1])
  
  return(name)
}


GenerateNetworkTitle <- function(param){
  
  title <- paste(c("Undirected","Directed")[param[["Directed"]]+1],
                 c("ident" = "", "abs" = "|weights|", "none" = "no weights",
                   "maxabs" = "max(|weights|)", "sum" = "sum(weights)",
                   "signedmaxabs" = "+/-max(|weights|)")[param[["Weights"]]],
                 c("unnormalized","row-normalized")[param[["Normalized"]]+1],
                 sep = ", ")
  
  return(title)
}


# Infer parameters used to generate network representation from file name

GetParamFromName <- function(name){
  
  param <- list("Directed" = !grepl("undirected", name),
                "Weights" = gsub("weights", "",
                                 strsplit(name, split = "_")[[1]][grep("weights",
                                                                       strsplit(name, split = "_")[[1]])]),
                "Normalized" = grepl("rownormalized", name))
  
  return(param)
  
}


# Read network files in external folders

ReadNetwork <- function(NET_FILE, NETLIBS){
  
  source(NETLIBS)
  
  network <- loadNetworkWithFilteringForSignificantPredictors(NET_FILE,
                                                              pValCutoff=-1,
                                                              localGeneCutoff = 50)
  net <- network$G
  rownames(net) <- network$responseGenes
  colnames(net) <- network$predictors
  rm(network); gc()
  net <- net[,-1]
  net <- net[rowSums(net) != 0,]
  net <- net[,colSums(net) != 0]
  
  return(net)
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


# Given a GO ID, returns the genes annotated to that term
GetGOGenes <- function(GOID){
  
  library(biomaRt)
  
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  genes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","go_id"),
                 filters = "go", values = GOID, mart = ensembl)
  
  return(genes)
}


# GOenrich: output from GetGOEnrich()
# col: color for bars
PlotGOEnrich <- function(GOenrich, col, title){
  
  library(ggplot2)
  
  p <- ggplot(GOenrich) +
    geom_bar(aes(x = Term, y = log2Enrichment, alpha = -log10(as.numeric(pval))),
             stat = "identity", fill = col) + guides(alpha = "none") +
    ggtitle(title) + coord_flip() + theme_classic() +
    theme(text = element_text(size = 20))
  
  return(p)
}



# Network functions -----------------------------------------------------------

# Convert adjacency matrix to igraph object

AdjToIgraph <- function(net, directed = T, weighted = T){
  
  library(igraph)
  
  if(directed){
    g <- igraph::graph_from_adjacency_matrix(t(net), mode = "directed",
                                             weighted = weighted)
  } else{
    stop("Not implemented for undirected")
  }
  
  return(g)
}


# Generate appropriate adjacency matrices

CreateNetworkMatrix <- function(param, net_graph, path){
  
  message("\nCreating network representation with parameters:")
  message(paste(names(param), param, collapse = "\n", sep = " = "))
  
  # Assumes param[["Directed]] == T
  if(param[["Weights"]] == "none"){
    message("Removing weights")
    net_graph <- remove.edge.attribute(net_graph, "weight")
  } else if(param[["Weights"]] == "abs"){
    message("Replacing weights by absolute value")
    E(net_graph)$weight <- abs(E(net_graph)$weight)
  }
  
  if(!param[["Directed"]]){
    message("Removing directionality")
    if(param[["Weights"]] == "none"){
      net_graph <- as.undirected(net_graph, mode = "collapse")
    } else{
      # need to handle weights
      net_graph <- as.undirected(net_graph, mode = "collapse", 
                                 edge.attr.comb = list(weight = switch(param[["Weights"]],
                                                                       "maxabs" = function(x) max(abs(x)),
                                                                       "signedmaxabs" = function(x) x[which.max(abs(x))],
                                                                       "sum" = function(x) sum(x))))
    }
  }
  
  # get adjacency matrix
  if(param[["Weights"]] == "none"){
    adj <- as_adjacency_matrix(net_graph, type = "both", sparse = FALSE)
  } else{
    adj <- as_adjacency_matrix(net_graph, type = "both", attr = "weight", sparse = FALSE)
  }
  
  if(param[["Normalized"]]){
    # row-normalize adjacency matrix
    message("Performing row-wise normalization by (weighted) degree")
    nf <- rowSums(abs(adj))
    nf[nf == 0] <- 1 # avoids dividing by 0, and keeps 0-only rows as they were
    adj <- sweep(adj, 1, nf, "/")
    rm(nf); gc()
  }
  
  # Save
  message("Saving network: ", GenerateNetworkName(param))
  saveRDS(adj, paste0(path,GenerateNetworkName(param),".rds"))
  
  return(NULL)
}


# Iterative Network Multiplication (like ADImpute but for all genes and not
# only dropouts)

IterativeNetMultiplication <- function(vector, network, max.iter = 50){
  
  # targets both in the network and vector
  t <- intersect(names(vector),rownames(network))
  # predictors that actually predict any of the targets AND are included in vector
  p <- intersect(names(vector), colnames(network)[Matrix::colSums(network[t, ]) != 0])
  
  network <- network[t,p]
  
  message("Starting network iterative imputation\n")
  
  i <- 1
  repeat {
    if ((max.iter != -1) & (i > max.iter)) {
      break
    }
    if ((i%%5) == 0) {
      message(paste("Iteration", i, "/", max.iter, "\n"))
    }
    
    new <- round(network %*% vector[p], 2)
    # expression = network coefficients * predictor expr.
    
    # Check convergence
    if (any(new[t,1] != vector[t])) {
      vector[t] <- new[t,1]
    } else {
      break
    }
    i <- i + 1
  }
  
  return(vector)
}



# Extend a non-square matrix to square (e.g. for a complete adjacency matrix
# that can be used to generate a graph)

MatrixToSquare <- function(mat){
  
  gg <- unique(c(colnames(mat),rownames(mat)))
  mat_ext <- matrix(0, nrow = length(gg), ncol = length(gg),
                    dimnames = list(gg, gg))
  mat_ext[rownames(mat),colnames(mat)] <- as.matrix(mat)
  
  return(mat_ext)
}



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
  
  message("Removing ",nrow(to_remove)," edges")
  
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
  
  p <- plot(graph, layout = switch(layout,
                                   "lgl" = layout_with_lgl,
                                   "fr" = layout_with_fr,
                                   "nicely" = layout_nicely),
       vertex.size = 5, edge.curved = .1, vertex.label.dist = 1,
       edge.arrow.size = .2, vertex.frame.color = NA,
       vertex.label.color = "grey50", vertex.label.family = "Helvetica",
       vertex.label.cex = .8, main = title, margin = 0)
  
  return(p)
}


SmoothNetworkRWR <- function(net, expr, alpha,
                             fakenet_file = "Outputs/Human_Network/stabsel/network_largest_cc_igraph.rds"){

  net_graph <- readRDS(fakenet_file)
  
  # Convert to "edge list"
  el <- as_edgelist(net_graph)
  rm(net_graph); gc()
  
  # Map expression values to network
  expr <- data.frame("gene_id" = rownames(expr), expr)
  netMap <- network_mapping(network = el, expr_mat = expr,
                            merge.by = "gene_id",
                            global = TRUE)
  
  # Replace altered network representation by original one
  netMap$G <- net[netMap$gene_names, netMap$gene_names]
  
  # smooth for all conditions
  smoothed <- network_smoothing(net = netMap$G,
                                mat_intensities = netMap$mat_intensities,
                                conditions = colnames(expr)[-1],
                                iter = 50, alpha = alpha)
  rownames(smoothed) <- netMap$gene_names
  
  return(smoothed)
}



# Plot ------------------------------------------------------------------------

# args should be a vector containing x variable, plot title and color

BarPlot <- function(plot.data, args){
  
  library(ggplot2)
  
  x <- args["x"]
  if("Title" %in% names(args)){
    title <-  args["Title"]
  } else{
    title <- NULL
  }
  if("Color" %in% names(args)){
    color <- args["Color"]
  }
  
  # keep general x axis limits
  xlims <- c(min(plot.data[[x]]), max(plot.data[[x]]))
  
  if(as.logical(args["Subset"]))
    plot.data <- subset(plot.data, subset = plot.data[[title]])
  
  p <- ggplot(plot.data) +
    geom_bar(aes(x = plot.data[[x]]), fill = color, color = color) +
    ggtitle(title) + xlab(x) +
    coord_cartesian(xlim = xlims) +
    theme_classic() + theme(text = element_text(size = 20))
  
  return(p)
}

DensityPlot <- function(plot.data, args){
  
  library(ggplot2)
  
  x <- args["x"]
  if("Title" %in% names(args)){
    title <-  args["Title"]
  } else{
    title <- NULL
  }
  if("Color" %in% names(args)){
    color <- args["Color"]
  } else{
    color <- "indianred"
  }
  
  # keep general x axis limits
  xlims <- c(min(plot.data[[x]]), max(plot.data[[x]]))
  
  if(as.logical(args["Subset"]))
    plot.data <- subset(plot.data, subset = plot.data[[title]])
  
  p <- ggplot(plot.data) +
    geom_density(aes(x = plot.data[[x]]), fill = color, color = color) +
    ggtitle(title) + xlab(x) +
    coord_cartesian(xlim = xlims) +
    theme_classic() + theme(text = element_text(size = 20),
                            axis.text.y = element_blank(),
                            axis.ticks.y = element_blank())
  
  return(p)
}


# Read / write ----------------------------------------------------------------

# Read file
ReadRDS <- function(file){
  message("Reading ",file)
  x <- readRDS(file)
  
  return(x)
}

# Write file
WriteRDS <- function(x, file){
  message("Writing ",file)
  saveRDS(x, file)
  
  return(NULL)
}


age_palette <- c("20-29" = "#126782", "30-39" = "#219EBC",
                 "50-59" = "#FD9E02", "60-69" = "#FB8500")
net_color <- "#1f73e0"
obs_color <- "olivedrab3"
trans_color <- "#06d6a0"
cis_color <- "#ff5a5f"
# tissue_palette <- c("Adipose-Subcutaneous" = "#77AD78",
#                     "Brain" = "#80A4ED",
#                     "Muscle-Skeletal" = "#EF767A")



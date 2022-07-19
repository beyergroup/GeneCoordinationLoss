net = "stabsel"
type = "network_largest_cc" # args[2]
weights = "none" # args[3]
matrix = "Adjacency" # args[4]
directed = T
normalize = F

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

if(directed){
  network <- readRDS("/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel/network_largest_cc_directed_adj_mat.rds")
} else{
  network <- readRDS("")
}


if(normalize){
  # correct network to remove topology bias effects - row-normalized adjacency
  nf <- rowSums(abs(network))
  nf[nf == 0] <- 1 # avoids dividing by 0, and keeps 0-only rows as they were
  network <- sweep(network, 1, nf, "/")
  rm(nf); gc()
}

tissue_files <- list.files("/data/public/adesous1/GeneCorrelation/Outputs/GTEx/logFC",
                           pattern = "sampled_meanLFC.rds", full.names = T)

smoothed <- list()
unsmoothed <- list() # just to compare to original LFCs

for(file in tissue_files){
  
  # read in tissue-specific GTEx subset for sample identification
  LFC <- readRDS(file)
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  smoothed[[tissue]] <- IterativeNetMultiplication(vector = LFC,
                                                   network = network)
  unsmoothed[[tissue]] <- LFC
  
  rm(LFC,tissue); gc()
}

smoothed <- do.call(cbind, smoothed)
unsmoothed <- do.call(cbind, unsmoothed)

# Get module membership -------------------------------------------------------

# read in decisions
decisions <- readRDS(paste0("Outputs/Human_Network/",net,"/Modules/",matrix,
                            "_weight",weights,"_height_eigen_decision.rds"))

# read in corresponding hc object
hc <- readRDS(paste0("Outputs/Human_Network/",net,"/Modules/",matrix,"_weight",
                     weights,"_complete_link_clustering_",
                     as.character(decisions$Eigenvectors),"_evectors.rds"))

# cut at chosen height
membership <- cutree(hc, h = decisions$Height)

# remove smaller than 10 genes
modules <- names(table(membership)[table(membership) >= 10])

rm(decisions,hc); gc()


# Summarize by module ---------------------------------------------------------

# Unsmoothed
activity <- data.frame(row.names = colnames(unsmoothed))
for(m in modules){
  genes <- names(which(membership == m))
  activity <- cbind.data.frame(activity,
                               as.data.frame(apply(unsmoothed[intersect(genes,rownames(unsmoothed)),],
                                                   2, mean)))
}
activity <- t(activity)
rownames(activity) <- modules

saveRDS(activity, paste0("Outputs/Human_Network/",net,
                         "/Modules/Module_Activity/NoNet_relative_activity_GTEx.rds"))

# Smoothed
activity <- data.frame(row.names = colnames(smoothed))
for(m in modules){
  genes <- names(which(membership == m))
  activity <- cbind.data.frame(activity,
                               as.data.frame(apply(smoothed[intersect(genes,rownames(smoothed)),],
                                                   2, mean)))
}
activity <- t(activity)
rownames(activity) <- modules

saveRDS(activity, paste0("Outputs/Human_Network/",net,
                         "/Modules/Module_Activity/",
                         c("weight","Adjacency_weight")[normalize+1],
                         weights,"_",c("undirected","directed")[directed+1],
                         "_mm50_relative_activity_GTEx.rds"))

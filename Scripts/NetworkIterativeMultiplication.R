# directed network, targets in rows and predictors in columns
network_directed <- readRDS("/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel/directed_adj_mat.rds")
# undirected version of the network (both edges take the max of the absolute values between the 2)
network_undirected <- readRDS("/data/public/adesous1/GeneCorrelation/Outputs/Human_Network/stabsel/undirected_maxabs_adj_mat.rds")

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


test_vector <- sample(c(0,1), size = 20000, replace = T)
names(test_vector) <- sample(rownames(network_directed), size = length(test_vector))
# first 10 entries not matched in network - function preserves them
names(test_vector)[1:10] <- paste0("Fake",1:10)

out1 <- IterativeNetMultiplication(vector = test_vector, network = network_directed)
out2 <- IterativeNetMultiplication(vector = test_vector, network = network_undirected)

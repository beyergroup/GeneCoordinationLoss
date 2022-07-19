
args = commandArgs(trailingOnly=TRUE)
net = args[1]
type = args[2]
weights = args[3]
matrix = args[4]

setwd("/data/public/adesous1/GeneCorrelation/")


# read in eigenvectors
ee <- readRDS(paste0("Outputs/Human_Network/",net,
                     "/Modules/",matrix,"_weight",weights,"_evectors_evalues_",type,".rds"))

N = c(
  ncol(ee$vectors),20000,
  15000,10000,5000,1000,500,100,50,20,10,5,
  NULL)

# base dist implementation is sooooo slow, use this instead:
euc_dist <- function(m) {mtm <- Matrix::tcrossprod(m); sq <- rowSums(m*m);  sqrt(outer(sq,sq,"+") - 2*mtm)}

for(n in N){
  
  message("Clustering on first ",n," eigenvectors")
  
  if(matrix == "Adjacency"){
    # take highest values from adjacency
    ee_ind <- 1:n
  } else if(matrix == "Laplacian"){
    # take lowest values from Laplacian
    ee_ind <- (1+(ncol(ee$vectors)-n)):ncol(ee$vectors)
  }
  d <- euc_dist(ee$vectors[,ee_ind])
  
  if(sum(is.na(d)) > 0){
    
    # if there are NAs, use integrated dist function to re-compute those distances
    na_ind <- which(is.na(d), arr.ind = TRUE)
    dists <- apply(na_ind, 1, function(x) dist(x = ee$vectors[c(x[1],x[2]),ee_ind]))
    for(i in 1:nrow(na_ind))
      d[na_ind[i,1],na_ind[i,2]] <- dists[i]
  }
  
  dd <- as.dist(d)
  
  message("Distances computed")
  
  for(met in c("complete","average")){
    hc <- hclust(dd, method = met)
    
    pdf(paste0("Plots/Human_Network/",net,"/Modules/",matrix,"_weight",weights,
               "_hclust_",n,"_",met,".pdf"), width = 20)
    plot(hc, labels = FALSE, main = c("stabsel" = "Indirect effects included",
                                      "stabsel_pcclasso" = "Indirect effects excluded")[net],
         sub = paste0(n," eigenvectors, ",met," linkage"))
    dev.off()
    saveRDS(hc, paste0("Outputs/Human_Network/",net,"/Modules/",matrix,"_weight",weights,
                       "_",met,"_link_clustering_",n,"_evectors.rds"))
    rm(hc); gc()
  }
  
  message("Clustering finished")
  
  rm(d,dd); gc()
}

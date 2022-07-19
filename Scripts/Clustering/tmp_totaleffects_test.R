# run pseudo-inversion on smaller subset of the network

# Moore-Penrose pseudo-inversion of (I-A)

args <- commandArgs(trailingOnly = TRUE)
matrix_file = args[1]
threshold = as.numeric(args[2])
out_dir = args[3]

message(threshold)

WDIR = "/data/public/adesous1/GeneCorrelation/"
setwd(WDIR)

.libPaths("Resources/Rlibs/R-3.4.4/")
library(MASS)

options(matprod="internal")

# Read in network matrix ------------------------------------------------------

message("Reading file ",matrix_file)
W <- readRDS(matrix_file)
message("Network matrix successfully read from file")
W <- W[1:100,1:100]


# Invert (I-A) ----------------------------------------------------------------

pinv <- ginv((diag(nrow(W)) - W), tol = threshold)
message("Pseudo-inverse computed")

# Save ------------------------------------------------------------------------

saveRDS(W %*% pinv,
        paste0(out_dir,gsub(".rds","_totaleffects.rds",
                            tail(strsplit(matrix_file,"/")[[1]],1))))
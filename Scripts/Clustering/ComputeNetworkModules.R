# Compute Modules (spectral clustering)

args <- commandArgs(trailingOnly = TRUE)
file = args[1]
INDIR = args[2]
OUTDIR = args[3]
PLOTDIR = args[4]

if(!(grepl("undirected",file))){
  stop("No directed representations accepted.")
}
if(grepl("signed",file)){
  stop("No negative weights accepted.")
}
if(!(grepl("largest_cc",file))){
  stop("Only largest connected component accepted.")
}

# Read in network
message("Reading in ",INDIR,file)
W <- readRDS(paste0(INDIR,file))

# Compute eigenvectors and eigenvalues of matrix
message("Computing eigenvectors and eigenvalues")
ee <- eigen(W, symmetric = TRUE)
rownames(ee$vectors) <- rownames(W)
message("Saving to ",OUTDIR)
saveRDS(ee, paste0(OUTDIR, gsub(".rds","_evectors_evalues.rds",file)))
message("Min value = ", min(abs(ee$values)))
sum(ee$values == min(abs(ee$values)))

# Plot
message("Plotting to ",PLOTDIR)
eigenvalues <- sort(ee$values, decreasing = FALSE)
pdf(paste0(PLOTDIR, gsub(".rds","_evalues.pdf",file)),
    height = 4)
plot(eigenvalues, main = gsub(".rds","",file))
dev.off()
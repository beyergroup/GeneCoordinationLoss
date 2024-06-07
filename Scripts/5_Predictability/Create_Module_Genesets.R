# Transform network modules to gene sets

NET = "stabsel_filtered_trans_largestCC"
OUT_DIR = "Outputs/5_Predictability/Age/Age_MaxSubset/SmoothedSlopes/GSEA/Inputs"
TISSUES = "main"

source("Scripts/functions.R")
source("Scripts/5_Predictability/params.R")

dir.create(OUT_DIR, recursive = T)

membership <- ReadRDS(paste0("Outputs/2_Clustering/",NET,"_",
                             "Adjacency_undirected","_",
                             "weightssum_rownormalized","_",
                             "reclustered_membership.rds"))

modules <- names(table(membership))[table(membership) >= 10]

sink(paste0(OUT_DIR,"/",NET,"_reclustered_genesets.gmt"))
for(module in modules){
  cat(paste(c(module,NA,names(which(membership == module))),
            collapse = "\t"))
  cat("\n")
}
sink()

# ----

# Format smoothed slopes into ranked lists

files <- list.files("Outputs/5_Predictability/Age/Age_MaxSubset/SmoothedSlopes",
                    pattern = "RWR_0.2")
for(file in files){
  smoothed <- ReadRDS(paste0("Outputs/5_Predictability/Age/Age_MaxSubset/SmoothedSlopes/",
                             file))
  # Write down the slopes
  write.table(smoothed[,"Slope"],
              file = paste0(OUT_DIR,"/",gsub(".rds","_slope.rnk",file)),
              sep = "\t", row.names = T, col.names = F, quote = F)
}


# Format original slopes into ranked lists

files <- list.files("Outputs/5_Predictability/Age/Age_MaxSubset",
                    pattern = "ageslope_well_predicted.rds", full.names = T)
slopes <- sapply(files, ReadRDS, simplify = F)
names(slopes) <- sapply(names(slopes),
                        function(x) strsplit(tail(strsplit(x, "/")[[1]],1),
                                             "_")[[1]][1])
rm(files); gc()

if(TISSUES != "all")
  slopes <- slopes[names(slopes) %in% MAIN_TISSUES]
common <- names(which(table(unlist(lapply(slopes,
                                          function(m) rownames(m)))) ==
                        length(slopes)))
slopes <- do.call(cbind, lapply(slopes, function(m) m[common,"Slope"]))

for(tissue in MAIN_TISSUES){
  # write.table(sort(slopes[,tissue], decreasing = T),
  #             file = paste0(OUT_DIR,"/",tissue,"_ageslope_original.rnk"),
  #             sep = "\t", row.names = T, col.names = F, quote = F)
  write.table(sort(slopes[[tissue]][,"Slope"], decreasing = T),
              file = paste0(OUT_DIR,"/",tissue,"_ageslope_all_original.rnk"),
              sep = "\t", row.names = T, col.names = F, quote = F)
}

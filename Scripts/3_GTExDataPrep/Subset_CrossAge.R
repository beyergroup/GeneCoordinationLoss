# Subset GTEx data

args <- commandArgs(trailingOnly = TRUE)
DATA_DIR = args[1]
SUBSET_DIR = args[2]
OUTDIR = args[3]

source("Scripts/functions.R")
source("Scripts/3_GTExDataPrep/params.R")

# metadata <- ReadRDS(paste0("Outputs/3_GTExDataPrep/metadata_restrictedaccess.rds"))
metadata <- ReadRDS(paste0(DATA_DIR,"/selected_annotation.rds"))

# subset metadata to fit sample filters
files <- list.files(DATA_DIR, pattern = "log2_kown_batch_corrected_gtex_")
samples <- sapply(paste0(DATA_DIR,"/",files),
                  function(f) rownames(ReadRDS(f)))
samples <- unlist(samples)
metadata <- subset(metadata, SAMPID %in% samples)
rm(samples); gc()


for(file in files[-grep("Brain|Kidney",files)]){
  
  tissue <- gsub(".rds","",tail(strsplit(file,"_")[[1]],1))
  
  # read in tissue subset
  gtex <- ReadRDS(paste0(DATA_DIR,"/",file))
  gtex <- t(gtex)
  
  # subset metadata to tissue and take into account samples that didnt't pass
  # filters
  m <- subset(metadata, SAMPID %in% colnames(gtex))
  
  # collect samples that were included in age groupings
  samples <- sapply(list.files(SUBSET_DIR, pattern = paste0(gsub(" ","",tissue),
                                                            ".+sampled_data"),
                               full.names = T),
                    function(f) colnames(ReadRDS(f)))
  samples <- unlist(samples)
  # exclude those samples from current subsetting
  m <- subset(m, !(SAMPID %in% samples))
  rm(samples); gc()
  
  if(all(table(m$AGE_GROUP)[c("30-39","40-49","50-59","60-69")] >= c(7,8,8,8))){
    samples <- mapply(function(age_group, n)
      sample(subset(m, AGE_GROUP == age_group)$SAMPID, size = n),
      c("30-39","40-49","50-59","60-69"),
      c(7,8,8,8))
    samples <- unlist(samples)
  } else{
    next
  }
  
  WriteRDS(gtex[,samples],
           paste0(OUTDIR,"/",gsub(" ","",tissue),"_sampled_data.rds"))
  
  rm(gtex,m,samples); gc()
}


# Brain
tissue <- "Brain"

gtex <- sapply(paste0(DATA_DIR,"/",files[grep(tissue, files, fixed = T)]),
               ReadRDS, simplify = F)
gtex <- do.call(rbind,gtex)
gtex <- t(gtex)

# subset metadata to tissue and take into account samples that didnt't pass
# filters
m <- subset(metadata, SAMPID %in% colnames(gtex))

# collect samples that were included in age groupings
samples <- sapply(list.files(SUBSET_DIR, pattern = paste0(gsub(" ","",tissue),
                                                          ".+sampled_data"),
                             full.names = T),
                  function(f) colnames(ReadRDS(f)))
samples <- unlist(samples)
# exclude those samples from current subsetting
m <- subset(m, !(SAMPID %in% samples))
rm(samples); gc()

samples <- mapply(function(age_group, n) sample(subset(m,
                                                       AGE_GROUP == age_group)$SAMPID,
                                                size = n),
                  c("30-39","40-49","50-59","60-69"),
                  c(7,8,8,8))
samples <- unlist(samples)

WriteRDS(gtex[,samples],
         paste0(OUTDIR,"/",gsub(" ","",tissue),"_sampled_data.rds"))

rm(gtex,m,samples); gc()

# Subset GTEx data

args <- commandArgs(trailingOnly = TRUE)
DATA_DIR = args[1]
SUBSET_DIR = args[2]

source("Scripts/functions.R")
source("Scripts/3_GTExDataPrep/params.R")

# metadata <- ReadRDS("Outputs/3_GTExDataPrep/metadata_restrictedaccess.rds")
metadata <- ReadRDS(paste0(DATA_DIR,"/selected_annotation.rds"))

# subset metadata to fit sample filters
files <- list.files(DATA_DIR, pattern = "log2_kown_batch_corrected_gtex_")
samples <- sapply(paste0(DATA_DIR,"/",files),
                  function(f) rownames(ReadRDS(f)))
samples <- unlist(samples)
metadata <- subset(metadata, SAMPID %in% samples)
rm(samples); gc()

# Join Brain subregions
metadata$Tissue <- metadata$SMTSD
metadata$Tissue[grepl("Brain",metadata$Tissue)] <- "Brain"

files <- list.files(DATA_DIR, pattern = "log2_kown_batch_corrected_gtex",
                    full.names = T)

for(tissue in MAIN_TISSUES){
  
  # read in tissue subset
  gtex <- sapply(files[grep(tissue, files, fixed = T)],
                 ReadRDS, simplify = F)
  gtex <- do.call(rbind,gtex)
  gtex <- t(gtex)
  
  # subset metadata to tissue and take into account samples that didnt't pass
  # filters
  m <- subset(metadata, SAMPID %in% colnames(gtex))
  
  # select number of samples to take from each age group (minimum number across
  # age groups excluding the 70s)
  t <- table(m$AGE_GROUP)
  n <- min(t[names(t)  != "70-79"])
  
  for(age_group in c("20-29","30-39","40-49","50-59","60-69","70-79")){
    
    if(tissue == "Brain"){
      
      s <- c()
      
      mm <- subset(m, AGE_GROUP == age_group)
      tt <- table(mm$SMTSD,mm$DONORID)
      tt <- tt[order(rowSums(tt)),order(colSums(tt))]
      
      for(stissue in rownames(tt)){
        
        if(grepl("Cortex|cortex|Cerebell|basal ganglia",stissue)){
          # pick 2 samples from these subtissues because they are similar to others
          d <- names(tt[stissue, tt[stissue,] != 0])
          d <- d[1:min(3,length(d))]
          s <- c(s, subset(mm, (DONORID %in% d) & (SMTSD == stissue))$SAMPID)
          tt <- tt[,c(setdiff(colnames(tt),d),d)]
        } else{
          d <- names(tt[stissue, tt[stissue,] != 0])
          d <- d[1:min(5,length(d))]
          s <- c(s, subset(mm, (DONORID %in% d) & (SMTSD == stissue))$SAMPID)
          tt <- tt[,c(setdiff(colnames(tt),d),d)]
        }
      }
    } else{
      s <- sample(subset(m, AGE_GROUP == age_group)$SAMPID,
                  size = min(nrow(subset(m, AGE_GROUP == age_group)),n))
    }
    
    # subset GTEx
    WriteRDS(gtex[,s], paste0(SUBSET_DIR,"/",
                              gsub(" ","",tissue),"_",age_group,
                              "_maxsampled_data.rds"))
  }
  rm(gtex,m,t,n); gc()
}




for(tissue in MAIN_TISSUES){
  
  # read in tissue subset
  gtex <- ReadRDS(paste0(DATA_DIR,"/log2_kown_batch_corrected_gtex_",
                         tissue,".rds"))
  gtex <- t(gtex)
  
  # subset metadata to tissue and take into account samples that didnt't pass
  # filters
  m <- subset(metadata, SAMPID %in% colnames(gtex))
  
  # select number of samples to take from each age group
  n <- min(sum(m$AGE < 40), sum(m$AGE >= 60))
  
  # young subset
  s <- sample(subset(m, AGE < 40)$SAMPID, size = n)
  WriteRDS(gtex[,s], paste0(SUBSET_DIR,"/",
                            gsub(" ","",tissue),"_young_maxsampled_data.rds"))
  
  # old subset
  s <- sample(subset(m, AGE >= 60)$SAMPID, size = n)
  WriteRDS(gtex[,s], paste0(SUBSET_DIR,"/",
                            gsub(" ","",tissue),"_old_maxsampled_data.rds"))
}

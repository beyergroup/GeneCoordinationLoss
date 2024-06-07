# Subset GTEx data

args <- commandArgs(trailingOnly = TRUE)
DATA_DIR = args[1]
SUBSET_DIR = args[2]
PLOT_DIR = args[3]

library(ggplot2)
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

# PCA on brain samples to determine which ones can be combined
gtex_brain <- sapply(paste0(DATA_DIR,"/",files[grep("Brain",files)]), ReadRDS)
gtex_brain <- do.call(rbind, gtex_brain)
pca <- prcomp(gtex_brain, center = T, scale. = T)
plot.data <- as.data.frame(pca$x[,1:2])
plot.data$SubTissue <- metadata[match(rownames(plot.data),metadata$SAMPID),"SMTSD"]
rm(gtex_brain,pca); gc()

ggplot(plot.data) +
  geom_point(aes(x = PC1, y = PC2, color = SubTissue)) +
  theme_classic() + theme(text = element_text(size = 20))
ggsave(PLOT_DIR,"/Brain_PCA.png")

# Join Brain subregions
metadata$Tissue <- metadata$SMTSD
metadata$Tissue[grepl("Brain",metadata$Tissue)] <- "Brain"

# Select age subset sizes
t <- table(metadata$Tissue, metadata$AGE_GROUP)
t <- t[rowSums(t[,1:5] >= SAMPLE_N) == 5,]


for(tissue in rownames(t)){
  
  # read in tissue subset
  gtex <- sapply(paste0(DATA_DIR,"/",files[grep(tissue, files, fixed = T)]),
                 ReadRDS, simplify = F)
  gtex <- do.call(rbind,gtex)
  gtex <- t(gtex)
  
  # subset metadata to tissue and take into account samples that didnt't pass
  # filters
  m <- subset(metadata, SAMPID %in% colnames(gtex))
  
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
          d <- d[1:min(2,length(d))]
          s <- c(s, subset(mm, (DONORID %in% d) & (SMTSD == stissue))$SAMPID)
          tt <- tt[,c(setdiff(colnames(tt),d),d)]
        } else{
          d <- names(tt[stissue, tt[stissue,] != 0])
          d <- d[1:min(3,length(d))]
          s <- c(s, subset(mm, (DONORID %in% d) & (SMTSD == stissue))$SAMPID)
          tt <- tt[,c(setdiff(colnames(tt),d),d)]
        }
      }
    } else{
      s <- sample(subset(m, AGE_GROUP == age_group)$SAMPID,
                  size = min(nrow(subset(m, AGE_GROUP == age_group)), SAMPLE_N))
    }
    
    if(length(s) < MIN_SAMPLES)
      next
    
    # subset GTEx
    WriteRDS(gtex[,s], paste0(SUBSET_DIR,"/",gsub(" ","",tissue),"_",
                              age_group,"_sampled_data.rds"))
  }
  rm(gtex,m,s); gc()
}


# for(tissue in MAIN_TISSUES){
#   
#   # read in tissue subset
#   gtex <- ReadRDS(paste0("/data/public/flopes/GTEx_AgeDeregulation_Paper/",
#                          "03-Known_Covariates_Batch_Correction/",
#                          "Outputs_Carolina_par/log2_kown_batch_corrected_gtex_",
#                          tissue,".rds"))
#   gtex <- t(gtex)
#   
#   # subset metadata to tissue and take into account samples that didnt't pass
#   # filters
#   m <- subset(metadata, SAMPID %in% colnames(gtex))
#   
#   # select number of samples to take from each age group
#   n <- min(sum(m$AGE < 40), sum(m$AGE >= 60))
#   
#   # young subset
#   s <- sample(subset(m, AGE < 40)$SAMPID, size = n)
#   WriteRDS(gtex[,s], paste0("Outputs/3_GTExDataPrep/Subset_Data/Max_Subset/",
#                             gsub(" ","",tissue),"_young_maxsampled_data.rds"))
#   
#   # old subset
#   s <- sample(subset(m, AGE >= 60)$SAMPID, size = n)
#   WriteRDS(gtex[,s], paste0("Outputs/3_GTExDataPrep/Subset_Data/Max_Subset/",
#                             gsub(" ","",tissue),"_old_maxsampled_data.rds"))
# }

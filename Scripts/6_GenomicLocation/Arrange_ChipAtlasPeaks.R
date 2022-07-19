# ChIP-seq peaks

AG = "LMNA"
TISSUE = "CDV"

celltypes <- list("Bld" = "Peripheral_blood_neutrophils",
                  "Epd" = "Human_epidermal_keratinocytes",
                  "CDV" = c("human_umbilical_vein_endothelial_cells",
                            "Human_Coronary_Artery_Endothelial_Cells_(HCAEC)",
                            "umbilical_vein_endothelial_cells"),
                  "Lng" = c("A549-CHIP-SMC1A","A549"))
# Digestive system only has a colorectal carcinoma cell line (SMC1A) or GP5d (SMC3)


peaks <- read.delim(paste0("Data/ChIP-Atlas/ChIPseq/Oth.",TISSUE,".50.",AG,
                           ".AllCell.bed"), header = F)
peaks <- peaks[-1,]
colnames(peaks)[1:3] <- c("Chromosome","StartPos","EndPos")
peaks$Chromosome <- gsub("chr","",peaks$Chromosome)

metadata <- sapply(peaks[,4], function(x) strsplit(x, ";")[[1]], simplify = F)
metadata <- lapply(metadata, function(x) gsub("%20","_",x))

# select right cell type
samples <- unlist(lapply(metadata,
                         function(m) {
                           if(any(grepl("source_name=",m))){
                            return(m[grep("source_name=",m)])
                          } else{
                            return(NA)
                          }}))
samples <- gsub("<br>source_name=","",samples)
names(samples) <- NULL

peaks <- peaks[samples %in% celltypes[[TISSUE]],]
metadata <- metadata[samples %in% celltypes[[TISSUE]]]

# experiment IDs (3 replicates)
table(unlist(lapply(metadata,
                    function(x){
                      if(length(intersect(grep("ID=",x),
                                          grep("SUBJECT_ID=",x,invert=T))) > 0){
                        return(x[intersect(grep("ID=",x),
                                           grep("SUBJECT_ID=",x,invert=T))])
                      } else{
                        return(NA)
                      }})))
peaks$Sample <- unlist(lapply(metadata,
                              function(x){
                                if(length(intersect(grep("ID=",x),
                                                    grep("SUBJECT_ID=",x,invert=T))) > 0){
                                  return(x[intersect(grep("ID=",x),
                                                     grep("SUBJECT_ID=",x,invert=T))])
                                } else{
                                  return(NA)
                                }}))
peaks$Tissue <- TISSUE
peaks$Target <- AG

saveRDS(peaks[,c("Chromosome","StartPos","EndPos","Sample","Tissue","Target")],
        paste0("Outputs/6_GenomicLocation/",AG,"_ChIP_",TISSUE,".rds"))

# Assign genes to TADs

files <- list.files("Data/hg38.TADs", pattern = "Schmitt2016")
files <- list.files("Data/hg38.TADs", pattern = "Leung_2015-raw_TADs.txt")


for(FILE in files){
  library(rCGH)
  source("Scripts/functions.R")
  
  
  chr.data <- rCGH::hg38
  rownames(chr.data)[c(23,24)] <- chr.data$chrom[c(23,24)] <- c("X","Y")
  
  all.genes <- ReadRDS("Outputs/6_GenomicLocation/gene_info.rds")
  
  # Read in TAD boundaries
  TADs <- read.delim(paste0("Data/hg38.TADs/",FILE), header = FALSE)
  colnames(TADs) <- c("Chromosome","StartPos","EndPos")
  TADs$Chromosome <- gsub("chr","",TADs$Chromosome)
  TADs$StartPosition <- TADs$StartPos + chr.data[TADs$Chromosome,"cumlen"]
  TADs$EndPosition <- TADs$EndPos + chr.data[TADs$Chromosome,"cumlen"]
  TADs$TAD_ID <- paste0("TAD",1:nrow(TADs))
  
  # Map genes to TADs
  all.genes$TAD_ID <- NA
  for(g in 1:nrow(all.genes)){
    if(any((all.genes$Position[g] >= TADs$StartPosition) &
           (all.genes$EndPosition[g] <= TADs$EndPosition))){
      all.genes$TAD_ID[g] <- TADs$TAD_ID[(all.genes$Position[g] >= TADs$StartPosition) &
                                           (all.genes$EndPosition[g] <= TADs$EndPosition)]
    }
  }
  
  TAD2gene <- merge(all.genes[,c("hgnc_symbol","Position","TAD_ID")],
                    TADs, by = "TAD_ID")
  WriteRDS(TAD2gene, paste0("Outputs/6_GenomicLocation/",
                            strsplit(FILE,"\\.")[[1]][1],"_genes.rds"))
}


library(ggplot2)


peaks <- read.delim("Data/ReMap/remap2022_peripheral-blood-neutrophil_nr_macs2_hg38_v1_0.bed",
                    header = F)

colnames(peaks)[1:4] <- c("Chromosome","StartPos","EndPos","Metadata")
peaks$Chromosome <- gsub("chr","",peaks$Chromosome)
peaks$Target <- sapply(peaks$Metadata, function(x) strsplit(x,":")[[1]][1])
peaks$MidPos <- rowMeans(cbind(peaks$StartPos,peaks$EndPos))


pdf("Plots/6_GenomicLocation/TADboundaries_neutrophil_ChIP.pdf",
    width = 20, height = 4)
for(chr in c(as.character(1:22),"X")){
  
  plot.data <- subset(peaks, Chromosome == chr, select = c("MidPos","Target"))
  
  print(ggplot(plot.data) +
          geom_point(aes(x = MidPos, y = Target), size = .2) +
          ggtitle(paste0("Chromosome ",chr)))
  rm(plot.data); gc()
}
dev.off()



WINDOW = 20

co_peaks <- data.frame("Chromosome" = c(), "MidPos" = c())

for(chr in c(as.character(1:22),"X")){
  
  message("Chromosome ",chr)
  
  # extract peaks of each target in the chromosome
  ctcf <- subset(peaks, (Chromosome == chr) & (Target == "CTCF"))$MidPos
  rad21 <- subset(peaks, (Chromosome == chr) & (Target == "RAD21"))$MidPos
  smc3 <- subset(peaks, (Chromosome == chr) & (Target == "SMC3"))$MidPos
  
  message(length(ctcf)," initial CTCF peaks")
  
  # for each CTCF peak check 50000 around for RAD21 and SMC3 peaks
  for(p in ctcf){
    if(any((rad21 > p-WINDOW) & (rad21 < p+WINDOW)) &&
       any((smc3 > p-WINDOW) & (smc3 < p+WINDOW))){
      co_peaks <- rbind.data.frame(co_peaks,
                                   data.frame("Chromosome" = chr,
                                              "MidPos" = p))
    }
  }
  
  message(nrow(subset(co_peaks, Chromosome == chr))," co-occuring peaks")
}


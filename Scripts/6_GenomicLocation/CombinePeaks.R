
TISSUE = "CDV"
TARGETS = c("CTCF","RAD21","SMC3","SMC1A")
WINDOW = 20

files <- list.files("Outputs/6_GenomicLocation", pattern = TISSUE,
                    full.names = T)
files <- files[grep(paste(TARGETS,collapse="|"),files)]

peaks <- sapply(files, readRDS, simplify = F)
peaks <- do.call(rbind, peaks)
peaks$MidPos <- rowMeans(peaks[,c("StartPos","EndPos")])


co_peaks <- data.frame("Chromosome" = c(), "MidPos" = c())

for(chr in c(as.character(1:22),"X")){
  
  message("Chromosome ",chr)
  
  message(nrow(subset(peaks, (Chromosome == chr) & (Target == "CTCF"))),
          " initial CTCF peaks")
  
  # iterate over middle of CTCF peaks
  for(p in subset(peaks, (Chromosome == chr) & (Target == "CTCF"))$MidPos){
    
    # check all other co-localizing proteins
    if(all(sapply(setdiff(unique(peaks$Target),"CTCF"),
                  function(co) any((subset(peaks, (Chromosome == chr) &
                                           (Target == co))$MidPos > p-WINDOW) &
                                   (subset(peaks, (Chromosome == chr) & 
                                           (Target == co))$MidPos < p+WINDOW))))){
      # all co-localizing proteins are present - save peak
      co_peaks <- rbind.data.frame(co_peaks,
                                   data.frame("Chromosome" = chr,
                                              "MidPos" = p))
    }
  }
  
  message(nrow(subset(co_peaks, Chromosome == chr))," co-occuring peaks")
}

saveRDS(co_peaks,paste0("Outputs/6_GenomicLocation/combined_peaks_",
                        TISSUE,".rds"))

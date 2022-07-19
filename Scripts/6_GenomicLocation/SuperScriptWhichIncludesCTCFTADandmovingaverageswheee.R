library(ggplot2)

tissues <- list("Artery-Tibial" = "Aorta", "Lung" = c("A549","Lung"))


# Chromosome data
chr.data <- rCGH::hg38
rownames(chr.data)[c(23,24)] <- chr.data$chrom[c(23,24)] <- c("X","Y")

# Well-predicted genes
well_predicted_genes <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")

# Gene data
all.genes <- ReadRDS("Outputs/6_GenomicLocation/gene_info.rds")

# Error age LFCs
cor_files <- list.files("Outputs/5_Predictability",
                        pattern = "_corLFCs_YvsO_adj_spearman.rds",
                        full.names = T)
cor_files <- cor_files[-grep("trans|cis",cor_files)]


LFC_THRE = .5
tissue = "Artery-Tibial"
t <- "CDV"

targets <- intersect(well_predicted_genes[[tissue]],
                     all.genes$hgnc_symbol)
data <- data.frame()
for(target in targets){
  data <- rbind.data.frame(data,
                           data.frame("Target" = target,
                                      "TargetPosition" = subset(all.genes, hgnc_symbol == target)$Position,
                                      "TargetRelPosition" = subset(all.genes, hgnc_symbol == target)$start_position,
                                      "Chromosome" = subset(all.genes, hgnc_symbol == target)$chromosome_name))
}

ageLFCs <- readRDS(cor_files[grepl(tissue,cor_files,fixed = T)])
data$TargetLFC <- ageLFCs[data$Target]
rm(ageLFCs)


# TAD boundaries in Lung

Lung_TADs <- read.delim("Data/hg38.TADs/Lung_Schmitt2016-raw_TADs.txt",
                        header = F)
colnames(Lung_TADs) <- c("Chromosome","StartPos","EndPos")
Lung_TADs$Chromosome <- gsub("chr","",Lung_TADs$Chromosome)
Lung_TADs$StartPosition <- Lung_TADs$StartPos + chr.data[Lung_TADs$Chromosome,"cumlen"]
Lung_TADs$EndPosition <- Lung_TADs$EndPos + chr.data[Lung_TADs$Chromosome,"cumlen"]

A549_TADs <- read.delim("Data/hg38.TADs/A549_raw-merged_TADs.txt", header = F)
colnames(A549_TADs) <- c("Chromosome","StartPos","EndPos")
A549_TADs$Chromosome <- gsub("chr","",A549_TADs$Chromosome)
A549_TADs$StartPosition <- A549_TADs$StartPos + chr.data[A549_TADs$Chromosome,"cumlen"]
A549_TADs$EndPosition <- A549_TADs$EndPos + chr.data[A549_TADs$Chromosome,"cumlen"]


# TAD boundaries in Artery
Artery_TADs <- read.delim("Data/hg38.TADs/Aorta_STL002_Leung_2015-raw_TADs.txt",
                          header = F)
colnames(Artery_TADs) <- c("Chromosome","StartPos","EndPos")
Artery_TADs$Chromosome <- gsub("chr","",Artery_TADs$Chromosome)
Artery_TADs$StartPosition <- Artery_TADs$StartPos + chr.data[Artery_TADs$Chromosome,"cumlen"]
Artery_TADs$EndPosition <- Artery_TADs$EndPos + chr.data[Artery_TADs$Chromosome,"cumlen"]


# CTCF peaks in Lung

# Lung_CTCF <- read.delim("Data/ReMap/remap2022_lung_nr_macs2_hg38_v1_0.bed",
#                         header = F)
# Lung_CTCF <- Lung_CTCF[,1:3]
# colnames(Lung_CTCF) <- c("Chromosome","StartPos","EndPos")
# Lung_CTCF$Chromosome <- gsub("chr","",Lung_CTCF$Chromosome)
# Lung_CTCF <- subset(Lung_CTCF, Chromosome %in% data$Chromosome)
# Lung_CTCF$MidPos <- rowMeans(cbind(Lung_CTCF$StartPos,Lung_CTCF$EndPos))
CTCF <- readRDS(paste0("Outputs/6_GenomicLocation/combined_peaks_",t,".rds"))


# Age-related LFC plot
data$Chromosome <- factor(as.character(data$Chromosome),
                          levels = c(as.character(1:22),"X","Y"))
CTCF$Chromosome <- factor(as.character(CTCF$Chromosome),
                          levels = c(as.character(1:22),"X"))
Lung_TADs$Chromosome <- factor(as.character(Lung_TADs$Chromosome),
                               levels = c(as.character(1:22),"X"))
A549_TADs$Chromosome <- factor(as.character(A549_TADs$Chromosome),
                               levels = c(as.character(1:22),"X"))
Artery_TADs$Chromosome <- factor(as.character(Artery_TADs$Chromosome),
                                 levels = c(as.character(1:22),"X"))

pdf(paste0("Plots/6_GenomicLocation/test_",tissue,"_copeaks.pdf"),
    width = 100, height = 5)
for(chr in c(as.character(1:22),"X")){
  print(ggplot(subset(data, (abs(TargetLFC) > LFC_THRE) & Chromosome == chr)) +
          geom_jitter(aes(x = TargetRelPosition, y = "Error LFCs",
                          color = TargetLFC),
                      height = 0.2, width = 0) +
          # geom_segment(data = subset(Lung_TADs, Chromosome == chr),
          #              aes(x = StartPos, y = "Lung TAD start",
          #                  xend = EndPos, yend = "Lung TAD end")) +
          # geom_segment(data = subset(A549_TADs, Chromosome == chr),
          #              aes(x = StartPos, y = "A549 TAD start",
          #                  xend = EndPos, yend = "A549 TAD end")) +
          geom_segment(data = subset(Artery_TADs, Chromosome == chr),
                       aes(x = StartPos, y = "Artery TAD start",
                           xend = EndPos, yend = "Artery TAD end")) +
          geom_jitter(data = subset(CTCF, Chromosome == chr),
                      aes(x = MidPos, y = "CTCF peaks"),
                      color = "red", shape = 4, size = .1,
                      height = 0.1, width = 0) +
          # facet_wrap(~ Chromosome, ncol = 1, scales = "free_x") +
          scale_color_gradient2(low = "#126782", mid = "grey",
                                high = "#FB8500",
                                limits = c(-max(abs(data$TargetLFC)),
                                           max(abs(data$TargetLFC)))) +
          scale_x_continuous(expand = c(0,0)) +
          xlab("Genomic position of target gene") +
          labs(color = "Error Age LFC") +
          ggtitle(paste0("|Error Age LFC| > ",LFC_THRE)) +
          theme_bw() + theme(axis.ticks = element_blank(),
                             text = element_text(size = 20),
                             legend.position = "bottom"))
}
dev.off()


# sliding window of LFCs

for(WINDOW in c(2,3,4)){
  
  data[,paste0("SlideW",WINDOW)] <- NA
  
  data <- data[order(data$TargetPosition),]
  for(i in (WINDOW+1):(nrow(data)-WINDOW)){
    
    # check that window is not crossing chromosome boundaries
    if(data[(i-WINDOW),"Chromosome"] == data[i,"Chromosome"] &
       (data[(WINDOW+i),"Chromosome"] == data[i,"Chromosome"])){
      
      data[i,paste0("SlideW",WINDOW)] <- mean(data[(i-WINDOW):(i+WINDOW),
                                                   "TargetLFC"])
    }
  }
}


# generate null distribution

shuffled <- replicate(1000, sample(data$TargetLFC))

for(WINDOW in c(1,2,3)){
  
  message(WINDOW)
  
  slides <- apply(shuffled, 2, function(r, data){
    
    slide <- rep(NA,WINDOW)
    for(i in (WINDOW+1):(nrow(data)-WINDOW)){
      if(data[(i-WINDOW),"Chromosome"] == data[i,"Chromosome"] &
         (data[(WINDOW+i),"Chromosome"] == data[i,"Chromosome"])){
        
        slide <- c(slide, mean(r[(i-WINDOW):(i+WINDOW)]))
      } else{
        slide <- c(slide, NA)
      }
    }
    slide <- c(slide,rep(NA,WINDOW))
    
    return(slide)
  }, data)
  
  avg <- rowMeans(slides)
  se <- apply(slides, 1, sd)
  
  data[,paste0("MeanRandW",WINDOW)] <- avg
  data[,paste0("BottomRandW",WINDOW)] <- avg-(1.96*se)
  data[,paste0("TopRandW",WINDOW)] <- avg+(1.96*se)
  
  rm(slides,avg,se); gc()
}


WINDOW = 3

pdf(paste0("Plots/6_GenomicLocation/test_",tissue,"_copeaks_movavg_W",WINDOW,".pdf"),
    width = 50, height = 5)
for(chr in c(as.character(1:22),"X")){
  
  annotation_y <- max(subset(data, Chromosome == chr,
                             select = paste0("SlideW",WINDOW)), na.rm = T)
  
  print(ggplot(subset(data, Chromosome == chr)) +
    # randomized data
    geom_line(aes_string(x = "TargetRelPosition",
                         y = paste0("MeanRandW",WINDOW)),
              color = "grey", linetype = "dashed") +
    geom_ribbon(aes_string(x = "TargetRelPosition",
                           ymin = paste0("BottomRandW",WINDOW),
                           ymax = paste0("TopRandW",WINDOW)),
                alpha = .1) +
    # moving average
    geom_line(aes_string(x = "TargetRelPosition",
                         y = paste0("SlideW",WINDOW))) +
    # # TADs from Lung tissue
    # geom_segment(data = subset(Lung_TADs, Chromosome == chr),
    #              aes(x = StartPos, xend = EndPos),
    #              y = annotation_y + (annotation_y/10),
    #              yend = annotation_y + (annotation_y/5)) +
    # # TADs from cell line matching CTCF data
    # geom_segment(data = subset(A549_TADs, Chromosome == chr),
    #              aes(x = StartPos, xend = EndPos),
    #              y = annotation_y + (3*annotation_y/10),
    #              yend = annotation_y + (2*annotation_y/5)) +
    geom_segment(data = subset(Artery_TADs, Chromosome == chr),
                 aes(x = StartPos, xend = EndPos),
                 y = annotation_y + (annotation_y/10),
                 yend = annotation_y + (annotation_y/5)) +
    # CTCF peaks
    geom_jitter(data = subset(CTCF, Chromosome == chr),
                aes(x = MidPos, y = annotation_y + (5*annotation_y/20)),
                color = "red", shape = 4, size = .1,
                height = (annotation_y/22), width = 0) +
    xlab("Genomic position (bp)") +
    ylab("Age-related error LFC") +
    ggtitle(paste("Chromosome",chr,tissue)) +
    scale_x_continuous(expand = c(0,0)) +
    # scale_y_continuous(expand = expansion(mult = c(0.05,0),
    #                                       add = c(0,(2*annotation_y/5)))) +
    scale_y_continuous(limits = c(NA, annotation_y + (6*annotation_y/20))) +
    theme_bw() + theme(text = element_text(size = 20)))
}
dev.off()

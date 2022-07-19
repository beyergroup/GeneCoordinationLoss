# Lung genes and TADs

library(rCGH)
library(biomaRt)
library(ggplot2)

TISSUE = "Lung"

chr.data <- rCGH::hg38
rownames(chr.data)[c(23,24)] <- chr.data$chrom[c(23,24)] <- c("X","Y")


# Well-predicted genes
well_predicted_genes <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")
global_targets <- unique(do.call(c,well_predicted_genes))

# Gene data
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
all.genes <- getBM(attributes = c("start_position","end_position",
                                  "hgnc_symbol","chromosome_name"),
                   filters = "hgnc_symbol", mart = mart,
                   values = global_targets)
all.genes <- all.genes[-grep("PATCH|CTG|MT",all.genes$chromosome_name),]
# remove genes with duplicate entries
all.genes <- all.genes[!(all.genes$hgnc_symbol %in% 
                           all.genes$hgnc_symbol[duplicated(all.genes$hgnc_symbol)]),]
all.genes$Position <- chr.data[all.genes$chromosome_name,"cumlen"] +
  all.genes$start_position
rownames(all.genes) <- all.genes$hgnc_symbol
rm(mart); gc()


# A549 cell line TADs
A549_TADs <- read.delim("Data/hg38.TADs/A549_raw-merged_TADs.txt", header = F)
colnames(A549_TADs) <- c("Chromosome","StartPos","EndPos")
A549_TADs$Chromosome <- gsub("chr","",A549_TADs$Chromosome)
A549_TADs$StartPosition <- A549_TADs$StartPos + chr.data[A549_TADs$Chromosome,"cumlen"]
A549_TADs$EndPosition <- A549_TADs$EndPos + chr.data[A549_TADs$Chromosome,"cumlen"]
A549_TADs$TAD_ID <- 1:nrow(A549_TADs)


# Map genes to TADs
tad_gene <- matrix(ncol = 2)

for(chr in c(as.character(1:22),"X")){
  
  tmp_tad <- subset(A549_TADs, Chromosome == chr)
  tmp_gen <- subset(all.genes, chromosome_name == chr)
  
  # iterate over TADs of chromosome
  for(i in 1:nrow(tmp_tad)){
    tmp_g <- tmp_gen$hgnc_symbol[(tmp_gen$start_position >= tmp_tad$StartPos[i]) &
                                   (tmp_gen$end_position <= tmp_tad$EndPos[i])]
    if(length(tmp_g) > 0){
      tad_gene <- rbind(tad_gene,
                        cbind("Genes" = tmp_g,
                              "TAD_ID" = tmp_tad$TAD_ID[i]))
    }
    rm(tmp_g); gc()
  }
  rm(tmp_tad,tmp_gen); gc()
}
tad_gene <- tad_gene[-1,]
any(duplicated(tad_gene[,1])) # should be F (each gene only belongs to 1 TAD)


# Age error LFCs
targets <- intersect(well_predicted_genes[[TISSUE]],
                     all.genes$hgnc_symbol)
ageLFCs <- readRDS(paste0("Outputs/5_Predictability/",TISSUE,
                          "_corLFCs_YvsO_adj_spearman.rds"))

data <- data.frame("Gene" = targets,
                   "AgeLFC" = ageLFCs[targets],
                   "Chromosome" = all.genes[targets,"chromosome_name"],
                   "Position" = all.genes[targets,"start_position"],
                   "TAD" = tad_gene[match(targets,tad_gene[,1]),2])
data <- data[!is.na(data$TAD),]
r <- c()
for(chr in unique(data$Chromosome)){
  d <- subset(data, Chromosome == chr)
  r <- c(r, order(d$Position))
}
data$Rank <- r
data$TAD <- factor(as.character(data$TAD),
                   levels = as.character(sort(as.numeric(unique(data$TAD)))))

pdf("Plots/6_GenomicLocation/TAD_boundaries_A549_Lung.pdf",
    height = 10, width = 20)
for(chr in c(as.character(1:22),"X")){
  print(ggplot(subset(data, Chromosome == chr)) +
      geom_point(aes(y = TAD, x = Rank, color = AgeLFC), size = .1) +
      scale_color_gradient2(low = "#126782", mid = "grey",
                            high = "#FB8500",
                            limits = c(-max(abs(data$AgeLFC)),
                                       max(abs(data$AgeLFC)))) +
      scale_x_continuous(expand = c(0,0)) +
      xlab("Gene order in chromosome") + ylab("TADs") +
      ggtitle(paste("Chromosome",chr)) +
      theme(text = element_text(size = 20),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()))
}
dev.off()

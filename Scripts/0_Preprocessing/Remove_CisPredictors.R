# Remove predictors in same chromosome arm

args <- commandArgs(trailingOnly = T)
NET = args[1]
WDIR = args[2]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(biomaRt)
library(ggplot2)
library(dplyr)
library(rCGH)
source("Scripts/functions.R")

# Load network ----------------------------------------------------------------
net <- ReadRDS(paste0("Outputs/",WDIR,"/",NET,"_network_Hs.rds"))
net <- as.matrix(net)
net <- MatrixToSquare(net)

# Gather genomic info data -----------------------------------------------------
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
all.genes <- getBM(attributes = c("start_position","end_position",
                                  "hgnc_symbol","chromosome_name"),
                   filters = "hgnc_symbol", mart = mart,
                   values = unique(c(rownames(net),colnames(net))))
all.genes <- all.genes[-grep("PATCH|CTG|MT",all.genes$chromosome_name),]
# remove genes with duplicate entries
all.genes <- all.genes[!(all.genes$hgnc_symbol %in% 
                           all.genes$hgnc_symbol[duplicated(all.genes$hgnc_symbol)]),]
# sum chromosome lengths to positions within chromosome
chr.data <- rCGH::hg38
chr.data$CentrPosition <- rowMeans(chr.data[,c("centromerStart","centromerEnd")])
rownames(chr.data)[c(23,24)] <- chr.data$chrom[c(23,24)] <- c("X","Y")
all.genes$Position <- all.genes$start_position
rownames(all.genes) <- all.genes$hgnc_symbol

# Check each predictor-target edge for membership to same chromosome arm -------
pairs <- data.frame()
for(target in rownames(net)){
  
  predictors <- intersect(names(net[target, net[target,] != 0]),
                          all.genes$hgnc_symbol)
  
  if(length(predictors) > 0){
    pairs <- rbind.data.frame(pairs,
                              data.frame("Target" = target,
                                         "Predictor" = predictors,
                                         "TargetChromosome" = all.genes[target,"chromosome_name"],
                                         "TargetPosition" = all.genes[target,"Position"],
                                         "PredictorChromosome" = all.genes[predictors,"chromosome_name"],
                                         "PredictorPosition" = all.genes[predictors,"Position"]))
  }
}
pairs$TargetArm <- paste0(pairs$TargetChromosome,"_",
                          as.character(as.numeric(pairs$TargetPosition >
                                                    chr.data[pairs$TargetChromosome,"CentrPosition"])+1))
pairs$PredictorArm <- paste0(pairs$PredictorChromosome,"_",
                             as.character(as.numeric(pairs$PredictorPosition >
                                                       chr.data[pairs$PredictorChromosome,"CentrPosition"])+1))

# Split network based on arm membership 
trans_net <- matrix(data = 0,
                  nrow = length(intersect(rownames(net),
                                          unique(pairs$Target))),
                  ncol = length(intersect(colnames(net),
                                          unique(pairs$Predictor))),
                  dimnames = list(intersect(rownames(net),
                                            unique(pairs$Target)),
                                  intersect(colnames(net),
                                            unique(pairs$Predictor))))

for(target in rownames(trans_net)){
  cur <- subset(pairs, Target == target)
  # cis predictors are genes in the same chromosome arm
  cis_pred <- subset(cur, PredictorArm == unique(TargetArm))$Predictor
  # trans predictors are genes whose position we know and are not cis
  trans_pred <- cur$Predictor[!(cur$Predictor %in% cis_pred)]
  trans_net[target,trans_pred] <- net[target,trans_pred]
}

WriteRDS(trans_net, paste0("Outputs/",WDIR,"/",NET,"_trans_network_Hs.rds"))
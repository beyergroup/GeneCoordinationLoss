# Distribution of genomic distance to the target gene

.libPaths("Resources/Rlibs/R-4.0.3/")
library(biomaRt)
library(ggplot2)
library(dplyr)
source("Scripts/functions.R")

# get list of well predicted genes
well_predicted_genes <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")
global_targets <- unique(do.call(c,well_predicted_genes))

# load network
net <- as.matrix(ReadRDS("Outputs/0_Preprocessing/stabsel_filtered_largestCC_network_Hs.rds"))
# reduce to well predicted targets (global)
net <- net[global_targets,]
net <- net[, colSums(net) != 0]
gc()

# get positions from biomaRt
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
all.genes <- getBM(attributes = c("start_position","end_position",
                                  "hgnc_symbol","chromosome_name"),
                   filters = "hgnc_symbol", mart = mart,
                   values = unique(c(global_targets,colnames(net))))
all.genes <- all.genes[-grep("PATCH|CTG|MT",all.genes$chromosome_name),]
# remove genes with duplicate entries
all.genes <- all.genes[!(all.genes$hgnc_symbol %in% 
                           all.genes$hgnc_symbol[duplicated(all.genes$hgnc_symbol)]),]
# sum chromosome lengths to positions within chromosome
chr_len <- c("1" = 247249719, "2" = 242951149, "3" = 199501827, "4" = 191273063,
             "5" = 180857866, "6" = 170899992, "7" = 158821424, "8" = 146274826,
             "9" = 140273252, "10" = 135374737, "11" = 134452384, "12" = 132349534,
             "13" = 114142980, "14" = 106368585, "15" = 100338915, "16" = 88827254,
             "17" = 78774742, "18" = 76117153, "19" = 63811651, "20" = 62435964,
             "21" = 46944323, "22" = 49691432, "X" = 154913754, "Y" = 57772954)
cum_chr_len <- c(0, sapply(2:24, function(i) sum(chr_len[1:(i-1)]),
                           simplify = T))
names(cum_chr_len) <- names(chr_len)
all.genes$Position <- cum_chr_len[all.genes$chromosome_name] +
  all.genes$start_position
rownames(all.genes) <- all.genes$hgnc_symbol


# Collect distances predictors-target
dists <- data.frame()
for(target in global_targets){
  
  predictors <- intersect(names(net[target, net[target,] != 0]),
                          all.genes$hgnc_symbol)
  
  if(length(predictors) > 0){
    dists <- rbind.data.frame(dists,
                              data.frame("Target" = target,
                                         "Predictor" = predictors,
                                         "TargetPosition" = all.genes[target,"Position"],
                                         "PredictorPosition" = all.genes[predictors,"Position"],
                                         "TargetChromosome" = all.genes[target,"chromosome_name"],
                                         "PredictorChromosome" = all.genes[predictors,"chromosome_name"]))
  }
}
dists$Distance <- abs(dists$TargetPosition - dists$PredictorPosition)


# Plot distribution of everyone's genomic position
ggplot(all.genes) +
  geom_density(aes(x = start_position)) +
  facet_wrap(~ chromosome_name, scales = "free")

# centromere is at the minimum gene density
centr <- c()
for(chr in c(as.character(1:22),"X")){
  d <- do.call(cbind,
               density(subset(all.genes,
                              chromosome_name == chr)$start_position)[c("x","y")])
  d <- d[floor(nrow(d)/5):ceiling(nrow(d)-(nrow(d)/5)),]
  centr <- c(centr, d[which.min(d[,"y"]),"x"])
}
names(centr) <- c(as.character(1:22),"X")

# Plot distribution of predictors-target distance (within chromosomes)
dists.chr <- data.frame()
for(target in unique(dists$Target)){
  cur <- subset(dists, Target == target)
  cur <- subset(cur, PredictorChromosome == unique(TargetChromosome))
  dists.chr <- rbind.data.frame(dists.chr,
                                data.frame("Distance" = cur$Distance,
                                           "Chromosome" = cur$TargetChromosome))
}

pdf("Plots/5_Predictability/pred_target_dists_by chr.pdf")
ggplot(dists.chr) +
  geom_density(aes(x = Distance)) +
  facet_wrap(~ Chromosome, scales = "free")
# plot(density(dists.chr),
#      main = "Genomic distance between target and predictor genes (in same chr.)",
#      xlab = "Genomic distance (bp)")
rm(cur, dists.chr); gc()


# Decide on a threshold
THRE = 20000000
abline(v = THRE)
dev.off()


# Plot genomic locations of predictors and targets using threshold

pdf("Plots/5_Predictability/genomic_positions_trans_cis.pdf",
    height = 20, width = 20)
ggplot(dists) +
  # trans
  geom_point(data = subset(dists, (PredictorChromosome != TargetChromosome) |
                             Distance >= THRE),
             aes(x = PredictorPosition, y = TargetPosition),
             size = .1, color = "#D36582") +
  # cis
  geom_point(data = subset(dists, (PredictorChromosome == TargetChromosome) &
                             Distance < THRE),
             aes(y = TargetPosition, x = PredictorPosition),
             size = .1, color = "#253C78") +
  geom_vline(xintercept = c(cum_chr_len,sum(chr_len)),
             linetype = "dashed") +
  geom_hline(yintercept = c(cum_chr_len,sum(chr_len)),
             linetype = "dashed") +
  scale_x_continuous(breaks = cumsum(chr_len) - (chr_len/2),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = cumsum(chr_len) - (chr_len/2),
                     expand = c(0,0)) +
  xlab("Genomic position of predictor gene(s)") +
  ylab("Genomic position of target gene") +
  theme_bw() + theme(axis.ticks = element_blank(),
                     text = element_text(size = 20))
dev.off()


# Apply threshold and split 
cis_net <- matrix(data = 0,
                  nrow = length(intersect(rownames(net),
                                          unique(dists$Target))),
                  ncol = length(intersect(colnames(net),
                                          unique(dists$Predictor))),
                  dimnames = list(intersect(rownames(net),
                                            unique(dists$Target)),
                                  intersect(colnames(net),
                                            unique(dists$Predictor))))
trans_net <- cis_net

for(target in rownames(cis_net)){
  cur <- subset(dists, Target == target)
  # cis predictors are genes in the same chromosome with dist below threshold
  cis_pred <- subset(cur, (PredictorChromosome == unique(TargetChromosome)) &
                       Distance < THRE)$Predictor
  # trans predictors are genes whose position we know and are not cis
  trans_pred <- cur$Predictor[!(cur$Predictor %in% cis_pred)]
  cis_net[target,cis_pred] <- net[target,cis_pred]
  trans_net[target,trans_pred] <- net[target,trans_pred]
}
intersect(which(cis_net != 0), which(trans_net != 0)) # no intersect, as supposed

# Save results
WriteRDS(cis_net, "Outputs/5_Predictability/cis_net.rds")
WriteRDS(trans_net, "Outputs/5_Predictability/trans_net.rds")

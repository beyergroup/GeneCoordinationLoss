# Distribution of genomic distance to the target gene

.libPaths("Resources/Rlibs/R-4.0.3/")
library(biomaRt)
library(ggplot2)
library(dplyr)
library(rCGH)
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
chr.data <- rCGH::hg38
chr.data$CentrPosition <- rowMeans(chr.data[,c("centromerStart","centromerEnd")])+chr.data$cumlen
rownames(chr.data)[c(23,24)] <- chr.data$chrom[c(23,24)] <- c("X","Y")
all.genes$Position <- chr.data[all.genes$chromosome_name,"cumlen"] +
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

# separate chromosomes into arms
dists$TargetArm <- paste0(dists$TargetChromosome,"_",
                          as.character(as.numeric(dists$TargetPosition >
                                                    chr.data[dists$TargetChromosome,"CentrPosition"])+1))
dists$PredictorArm <- paste0(dists$PredictorChromosome,"_",
                             as.character(as.numeric(dists$PredictorPosition >
                                                       chr.data[dists$PredictorChromosome,"CentrPosition"])+1))

# Plot genomic locations of predictors and targets in same chromosome arm or not
pdf("Plots/5_Predictability/genomic_positions_trans_cis_chrarm.pdf",
    height = 20, width = 20)
ggplot(dists) +
  # cis
  geom_point(data = subset(dists, TargetArm == PredictorArm),
             aes(x = PredictorPosition, y = TargetPosition),
             size = .1, color = cis_color) +
  # trans
  geom_point(data = subset(dists, TargetArm != PredictorArm),
             aes(x = PredictorPosition, y = TargetPosition),
             size = .1, color = trans_color) +
  geom_vline(xintercept = chr.data$cumlen, linetype = "dashed") +
  geom_hline(yintercept = chr.data$cumlen, linetype = "dashed") +
  scale_x_continuous(breaks = chr.data$cumlen + (chr.data$len/2),
                     expand = c(0,0), labels = chr.data$chrom) +
  scale_y_continuous(breaks = chr.data$cumlen + (chr.data$len/2),
                     expand = c(0,0), labels = chr.data$chrom) +
  xlab("Genomic position of predictor gene(s)") +
  ylab("Genomic position of target gene") +
  theme_bw() + theme(axis.ticks = element_blank(),
                     text = element_text(size = 20))
dev.off()



# Split network based on arm membership 
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
  # cis predictors are genes in the same chromosome arm
  cis_pred <- subset(cur, PredictorArm == unique(TargetArm))$Predictor
  # trans predictors are genes whose position we know and are not cis
  trans_pred <- cur$Predictor[!(cur$Predictor %in% cis_pred)]
  cis_net[target,cis_pred] <- net[target,cis_pred]
  trans_net[target,trans_pred] <- net[target,trans_pred]
}
intersect(which(cis_net != 0), which(trans_net != 0)) # no intersect, as supposed

# Save results
WriteRDS(cis_net, "Outputs/5_Predictability/cis_net_chrarm.rds")
WriteRDS(trans_net, "Outputs/5_Predictability/trans_net_chrarm.rds")


# Plot number of edges and total weights per target
sum(cis_net != 0)
sum(trans_net != 0)

pdf("Plots/5_Predictability/cis_trans_edgeweights.pdf")
plot(density(cis_net[cis_net != 0]), main = "Weight distribution - cis edges")
plot(density(trans_net[trans_net != 0]), main = "Weight distribution- trans edges")
dev.off()

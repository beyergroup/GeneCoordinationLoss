# Gather information on genomic location of genes and their predictors
# Join info on age-related changes in expression and predictability

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

# get error age LFCs
cor_files <- list.files("Outputs/5_Predictability",
                        pattern = "_corLFCs_YvsO.rds",
                        full.names = T)

# get expression age LFCs
expression_files <- list.files("Outputs/GTEx/AgeDE",
                               pattern = "_ageDE.rds",
                               full.names = T)


plot.data <- data.frame()

for(tissue in names(well_predicted_genes[-4])){
  
  targets <- intersect(well_predicted_genes[[tissue]],
                       all.genes$hgnc_symbol)
  
  data <- data.frame()
  
  for(target in targets){
    
    predictors <- intersect(names(net[target, net[target,] != 0]),
                            all.genes$hgnc_symbol)
    
    if(length(predictors) > 0){
      data <- rbind.data.frame(data,
                               data.frame(
                                 "Target" = target,
                                 "Predictor" = predictors,
                                 "TargetPosition" = subset(all.genes, hgnc_symbol == target)$Position,
                                 "PredictorPosition" = all.genes[predictors,"Position"],
                                 "TargetChromosome" = subset(all.genes, hgnc_symbol == target)$chromosome_name,
                                 "PredictorChromosome" = all.genes[predictors,"chromosome_name"]))
    }
  }
  
  ageLFCs <- readRDS(cor_files[grepl(tissue,cor_files,fixed = T)])
  
  data$TargetLFC <- ageLFCs[data$Target]
  data$PredictorLFC <- ageLFCs[data$Predictor]
  
  exprLFCs <- readRDS(expression_files[grepl(tissue,cor_files,fixed = T)])
  exprLFCs <- limma::topTable(fit = exprLFCs, coef = "AgeGroupOld",
                              number = nrow(exprLFCs))
  exprLFCs <- DataENSGToSymbol(exprLFCs, remove_dup = T)
  data$TargetExprLFC <- exprLFCs[data$Target,"logFC"]
  data$TargetExprFDR <- exprLFCs[data$Target,"adj.P.Val"]
  
  data$Tissue <- tissue
  plot.data <- rbind.data.frame(plot.data, data)
  rm(data); gc()
  
}

WriteRDS(plot.data, "Outputs/5_Predictability/genomic_location_info.rds")

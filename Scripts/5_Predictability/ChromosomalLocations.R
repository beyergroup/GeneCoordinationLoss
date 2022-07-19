# get genomic location table

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
                        pattern = "_corLFCs_YvsO_adj_spearman.rds",
                        full.names = T)
cor_files <- cor_files[-grep("trans|cis",cor_files)]

# get expression age LFCs
expression_files <- list.files("Outputs/GTEx/AgeDE",
                               pattern = "_ageDE.rds",
                               full.names = T)

LFC_THRE = .5



pdf("Plots/5_Predictability/genomic_positions_adj_spearman_05thre.pdf",
    height = 22, width = 20)
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
  
  cur.plot.data <- subset(data, abs(TargetLFC) > LFC_THRE,
                          select = c("Target","TargetPosition",
                                     "TargetLFC","TargetExprLFC",
                                     "TargetExprFDR","TargetChromosome"))
  cur.plot.data <- distinct(cur.plot.data)
  cur.plot.data$Tissue <- tissue
  plot.data <- rbind.data.frame(plot.data, cur.plot.data)
  rm(cur.plot.data); gc()
  
  print(ggplot(data) +
          geom_point(data = subset(data, abs(data$TargetLFC) <= LFC_THRE),
                     aes(y = TargetPosition, x = PredictorPosition,
                         color = TargetLFC), size = .1) +
          scale_color_gradient2(low = "#126782", mid = "grey",
                                high = "#FB8500",
                                limits = c(-max(abs(ageLFCs)),
                                           max(abs(ageLFCs)))) +
          geom_point(data = subset(data, abs(data$TargetLFC) > LFC_THRE),
                     aes(x = PredictorPosition, y = TargetPosition,
                         color = TargetLFC)) +
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
          labs(color = "Age LFC\n(target gene)") +
          ggtitle(tissue) +
          theme_bw() + theme(axis.ticks = element_blank(),
                             text = element_text(size = 20),
                             legend.position = "bottom"))
  
  print(ggplot(subset(data, abs(data$TargetLFC) > LFC_THRE)) +
          scale_color_gradient2(low = "#126782", mid = "grey",
                                high = "#FB8500",
                                limits = c(-max(abs(ageLFCs)),
                                           max(abs(ageLFCs)))) +
          geom_point(data = subset(data, abs(data$TargetLFC) > LFC_THRE),
                     aes(x = PredictorPosition, y = TargetPosition,
                         color = TargetLFC)) +
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
          labs(color = "Age LFC\n(target gene)") +
          ggtitle(gsub("(","\n(",tissue,fixed=T)) +
          theme_bw() + theme(axis.ticks = element_blank(),
                             text = element_text(size = 20),
                             legend.position = "bottom"))
  
  print(ggplot(subset(data, abs(data$TargetLFC) > LFC_THRE)) +
          scale_color_gradient2(low = "#126782", mid = "grey",
                                high = "#FB8500",
                                limits = c(-max(abs(ageLFCs)),
                                           max(abs(ageLFCs)))) +
          geom_point(data = subset(data, abs(data$TargetLFC) > LFC_THRE),
                     aes(x = PredictorPosition, y = TargetPosition,
                         color = PredictorLFC)) +
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
          labs(color = "Age LFC\n(predictor gene)") +
          ggtitle(gsub("(","\n(",tissue,fixed=T)) +
          theme_bw() + theme(axis.ticks = element_blank(),
                             text = element_text(size = 20),
                             legend.position = "bottom"))
}
dev.off()


# plot.data <- subset(plot.data, abs(TargetLFC) < 5)

pdf("Plots/5_Predictability/genomic_positions_alltissues_adj_spearman_05thre.pdf",
    height = 8, width = 20)

ggplot(plot.data) +
  geom_jitter(aes(x = TargetPosition, y = Tissue, color = TargetLFC),
              height = 0.2, width = 0) +
  geom_vline(xintercept = c(cum_chr_len,sum(chr_len)),
             linetype = "dashed") +
  scale_color_gradient2(low = "#126782", mid = "grey",
                        high = "#FB8500",
                        limits = c(-max(abs(ageLFCs)),
                                   max(abs(ageLFCs)))) +
  scale_x_continuous(breaks = cumsum(chr_len) - (chr_len/2),
                     expand = c(0,0)) +
  xlab("Genomic position of target gene") +
  labs(color = "Error Age LFC") + ggtitle(paste0("|Error Age LFC| > ",LFC_THRE)) +
  theme_bw() + theme(axis.ticks = element_blank(),
                     text = element_text(size = 20),
                     legend.position = "bottom")

ggplot(subset(plot.data, (TargetExprFDR < 0.05) & (TargetExprLFC > 0))) +
  geom_jitter(aes(x = TargetPosition, y = Tissue, color = TargetLFC),
              height = 0.2, width = 0) +
  geom_vline(xintercept = c(cum_chr_len,sum(chr_len)),
             linetype = "dashed") +
  scale_color_gradient2(low = "#126782", mid = "grey",
                        high = "#FB8500",
                        limits = c(-max(abs(ageLFCs)),
                                   max(abs(ageLFCs)))) +
  scale_x_continuous(breaks = cumsum(chr_len) - (chr_len/2),
                     expand = c(0,0)) +
  xlab("Genomic position of target gene") +
  labs(color = "Age LFC") +
  ggtitle(paste0("|Error Age LFC| > ",LFC_THRE,", up-regulated w/ age")) +
  theme_bw() + theme(axis.ticks = element_blank(),
                     text = element_text(size = 20),
                     legend.position = "bottom")

ggplot(subset(plot.data, (TargetExprFDR < 0.05) & (TargetExprLFC < 0))) +
  geom_jitter(aes(x = TargetPosition, y = Tissue, color = TargetLFC),
              height = 0.2, width = 0) +
  geom_vline(xintercept = c(cum_chr_len,sum(chr_len)),
             linetype = "dashed") +
  scale_color_gradient2(low = "#126782", mid = "grey",
                        high = "#FB8500",
                        limits = c(-max(abs(ageLFCs)),
                                   max(abs(ageLFCs)))) +
  scale_x_continuous(breaks = cumsum(chr_len) - (chr_len/2),
                     expand = c(0,0)) +
  xlab("Genomic position of target gene") +
  labs(color = "Age LFC") +
  ggtitle(paste0("|Error Age LFC| > ",LFC_THRE,", down-regulated w/ age")) +
  theme_bw() + theme(axis.ticks = element_blank(),
                     text = element_text(size = 20),
                     legend.position = "bottom")

dev.off()


# kNN-like assignment of LFCs as positive or negative
contingency <- list("4" = list(), "6" = list(), "8" = list())
pvalues <- matrix(nrow = 3, ncol = length(unique(plot.data$Tissue)),
                  dimnames = list(c("4","6","8"),unique(plot.data$Tissue)))
for(tissue in unique(plot.data$Tissue)){
  
  cur.plot.data <- subset(plot.data, Tissue == tissue)
  
  for(k in as.numeric(names(contingency))){
    
    predictions <- c()
    
    for(chr in unique(cur.plot.data$TargetChromosome)){
      
      tmp <- subset(cur.plot.data, TargetChromosome == chr)
      tmp <- tmp[order(tmp$TargetPosition),]
      
      if(nrow(tmp) < k+1)
        next
      
      p <- c()
      
      for(i in ((k/2)+1):(nrow(tmp)-(k/2))){
        p <- c(p, c("Neg","Pos")[1+as.numeric(median(tmp$TargetLFC[(i-(k/2)):(i+(k/2))] > 0))])
      }
      names(p) <- tmp$Target[((k/2)+1):(nrow(tmp)-(k/2))]
      predictions <- c(predictions,p)
      
      rm(tmp,p); gc()
    }
    
    observations <- c("Neg","Pos")[1+as.numeric(subset(cur.plot.data, Target %in% names(predictions))$TargetLFC > 0)]
    contingency[[as.character(k)]][[tissue]] <- table(predictions,observations)
    
    if((length(unique(predictions)) > 1) & (length(unique(observations)) > 1)){
      ftest <- fisher.test(x = as.factor(predictions),y = as.factor(observations))
      pvalues[as.character(k),tissue] <- ftest$p.value
    }
  }
}



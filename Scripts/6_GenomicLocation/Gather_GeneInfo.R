# Gather gene annotation data

library(rCGH)
library(biomaRt)
source("Scripts/functions.R")

chr.data <- rCGH::hg38
rownames(chr.data)[c(23,24)] <- chr.data$chrom[c(23,24)] <- c("X","Y")

# Well-predicted genes
well_predicted_genes <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")
global_targets <- unique(do.call(c,well_predicted_genes))
rm(well_predicted_genes); gc()

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
all.genes$EndPosition <- chr.data[all.genes$chromosome_name,"cumlen"] +
  all.genes$end_position
rownames(all.genes) <- all.genes$hgnc_symbol
rm(mart); gc()

WriteRDS(all.genes, "Outputs/6_GenomicLocation/gene_info.rds")

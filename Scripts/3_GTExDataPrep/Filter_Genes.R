# Filter out genes

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR = args[1]

librray(ggplot2)
source("Scripts/functions.R")

#Get GTEx original rownames(genenames):
GTEx_rownames <- rownames(readRDS(paste0(OUT_DIR,"/sampleFiltered_gtex.rds")));gc()
ensembl.genes <- cleanid(GTEx_rownames) #Remove the .version on genes ensemble id.
rm(GTEx_rownames);gc()
#Convert from ensembl.gene to gene.symbol
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("GENEID","SYMBOL","GENEBIOTYPE","SEQNAME"))
SaveRDS(geneIDs, file = paste0(OUT_DIR,"/original_geneIDs_GTEx_rownames.rds"))


####Select Genes By Biotype--------------------------------------------------------------------
#Check frequency table of gene's biotype:
GeneBiotypeCount<-as.data.frame(table(geneIDs[,3], dnn = "GENEBIOTYPE"))
GeneBiotypeCount<-GeneBiotypeCount[order(GeneBiotypeCount$Freq),]
GeneBiotypeCount<-GeneBiotypeCount[nrow(GeneBiotypeCount):1,]
rownames(GeneBiotypeCount) = NULL
SaveRDS(GeneBiotypeCount, file = paste0(OUT_DIR,"/original_GeneBiotypeCount.rds"))
rm(GeneBiotypeCount);gc()

#Check frequency table of gene's chromossome:
ChromossomeCount<-as.data.frame(table(geneIDs[,4], dnn = "Chromossome"))
ChromossomeCount<-ChromossomeCount[order(ChromossomeCount$Freq),]
ChromossomeCount<-ChromossomeCount[nrow(ChromossomeCount):1,]
rownames(ChromossomeCount) = NULL
SaveRDS(ChromossomeCount, file = paste0(OUT_DIR,"/original_ChromossomeCount.rds"))
rm(ChromossomeCount);gc()

#Select Genes By Biotype
geneIDs <- readRDS(paste0(OUT_DIR,"/original_geneIDs_GTEx_rownames.rds"))
SelectedGenesByBiotype <- geneIDs[
  geneIDs[,3]=="protein_coding"|
    geneIDs[,3]=="lincRNA"|
    geneIDs[,3]=="snRNA"|
    geneIDs[,3]=="miRNA"|
    geneIDs[,3]=="snoRNA",]
rownames(SelectedGenesByBiotype) = NULL
SaveRDS(SelectedGenesByBiotype, file = paste0(OUT_DIR,"/01-SelectedGenesByBiotype.rds"))
sampleFiltered_gtex <- readRDS(paste0(OUT_DIR,"/sampleFiltered_gtex.rds"))
sampleFiltered_gtex <- sampleFiltered_gtex[ensembl.genes %in% SelectedGenesByBiotype$GENEID,]
SaveRDS(sampleFiltered_gtex, paste0(OUT_DIR,"/sample_and_gene_biotype_filtered_gtex.rds"))
rm(sampleFiltered_gtex);gc()
#----------------------------------------------------------------------------------------------
  

####Select Genes By Value--------------------------------------------------------------------
#Check Top 50 gene's average expression in a density plot:

gtex <- readRDS(paste0(OUT_DIR,"/sample_and_gene_biotype_filtered_gtex.rds"))
means <- data.frame("Means" = apply(gtex, 1, function(x) mean(tail(sort(x),50))))
THRESHOLD_MEANS <- 100
sp <- ggplot(means, aes(x = Means)) +
  geom_density(color = "darkblue", fill = "lightblue") +
  geom_vline(xintercept = THRESHOLD_MEANS) +
  labs(x = "Genes Mean Expression") +
  scale_x_continuous(trans = 'log10') +
  labs(x = "Genes LOG10 Mean DESeq2 Normalized Expression")
ggsave(paste0(PLOT_DIR,"/means_log_prefiltering.png"),plot = sp)
gc()
SelectedGenesByMean <- as.data.frame(means[means[,"Means"]>=100,])
colnames(SelectedGenesByMean) <- colnames(means)
rownames(SelectedGenesByMean) = rownames(means)[means[,"Means"]>=100]
SaveRDS(SelectedGenesByMean, file = paste0(OUT_DIR,"/02-SelectedGenesByMean.rds"))
gtex <- gtex[rownames(gtex) %in% rownames(SelectedGenesByMean),]
SaveRDS(gtex, paste0(OUT_DIR,"/sample_and_gene_biotype_and_mean_filtered_gtex.rds"))
rm(gtex);gc()
#----------------------------------------------------------------------------------------------


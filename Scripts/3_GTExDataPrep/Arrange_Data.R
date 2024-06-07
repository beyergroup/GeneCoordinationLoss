
args <- commandArgs(trailingOnly = TRUE)
DATA_DIR = args[1]
OUT_DIR = args[2]

library(CePa)

## 00 - Data loading

sample_annotations <- read.delim(paste0(DATA_DIR,"/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"))
donor_annotations <- read.delim(paste0(DATA_DIR,"/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"))

# get donor ID from sample ID
sample_annotations$SUBJID <- sapply(sample_annotations$SAMPID,
  function(x) paste(strsplit(x, "-")[[1]][1:2], collapse = "-"))
annotation <- merge(sample_annotations, donor_annotations)
rm(sample_annotations,donor_annotations); gc()
SaveRDS(annotation, paste0(OUT_DIR,"/sample_donor_merged_annotation.rds"))


gtex <- read.gct(paste0(DATA_DIR,"/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"))

# reduce annotation to GTEx samples and apply the same order
colnames(gtex) <- gsub(x = colnames(gtex), pattern = "\\.", replacement = "-")
annotation <- annotation[match(colnames(gtex),annotation$SAMPID),]
rownames(annotation) <- NULL

saveRDS(gtex, paste0(OUT_DIR,"/raw_gtex.rds"))
saveRDS(annotation, paste0(OUT_DIR,"/annotation.rds"))


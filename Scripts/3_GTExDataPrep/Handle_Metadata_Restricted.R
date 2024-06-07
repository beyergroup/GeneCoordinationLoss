# Restricted access metadata

subj <- read.delim("/cellnet/AgeingClocks/Data/raw/gtex/GTEx/72188/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/PhenotypeFiles/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt",
                   sep = "\t", skip = 10)
samp <- read.delim("/cellnet/AgeingClocks/Data/raw/gtex/GTEx/72188/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/PhenotypeFiles/phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Attributes.GRU.txt",
                   sep = "\t", skip = 10)

samp$DONORID <- sapply(samp$SAMPID, function(x) paste(strsplit(x, "-")[[1]][1:2],
                                                      collapse = "-"))
all(samp$DONORID %in% subj$SUBJID)

# merge samples and donors metadata
metadata <- merge(x = samp, y = subj, by.x = "DONORID", by.y = "SUBJID")
rm(samp,subj); gc()

# translate sex info into readable
metadata$SEX <- c("Male","Female")[as.numeric(metadata$SEX)]

# create age groups from numeric ages
metadata$AGE_GROUP <- "20-29"
metadata$AGE_GROUP[as.numeric(metadata$AGE) > 29] <- "30-39"
metadata$AGE_GROUP[as.numeric(metadata$AGE) > 39] <- "40-49"
metadata$AGE_GROUP[as.numeric(metadata$AGE) > 49] <- "50-59"
metadata$AGE_GROUP[as.numeric(metadata$AGE) > 59] <- "60-69"
metadata$AGE_GROUP[as.numeric(metadata$AGE) > 69] <- "70-79"

# restrict only to RNA-seq data
metadata <- subset(metadata, SMAFRZE == "RNASEQ")

saveRDS(metadata, "Outputs/3_GTExDataPrep/metadata_restrictedaccess.rds")


table(metadata$SMTSD, metadata$AGE_GROUP)

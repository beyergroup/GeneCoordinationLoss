# Split by gender

metadata <- readRDS("Outputs/3_GTExDataPrep/metadata_restrictedaccess.rds")

metadata <- metadata[-grep("Cells",metadata$SMTSD),]
metadata <- droplevels.data.frame(metadata)

gtex <- readRDS("/cellnet/GeneCorrelation/Francisco/GTEx_analysis/04-Known_Covariates_Batch_Correction/Outputs/log2_batch_corrected_gtex.rds")
gtex <- t(gtex)

metadata <- subset(metadata, SAMPID %in% colnames(gtex))

# 150 samples per gender
t <- table(metadata$SMTSD, metadata$SEX)
gender_tissues <- rownames(t)[rowSums(t > 150) == ncol(t)]
rm(t); gc()

# 40 samples per age group
m <- subset(metadata, AGE_GROUP != "70-79")
t <- table(m$SMTSD, m$AGE_GROUP)
age_tissues <- rownames(t)[rowSums(t > 40) == ncol(t)]
rm(m,t); gc()


# smokers vs non-smokers
m <- subset(metadata, SMTSD == "Lung")
table(m$MHSMKSTS, m$AGE_GROUP)
# can set 20 samples per age group and "smoker"/"non-smoker" category

# donor in a ventilator
table(m$DTHVNT, m$AGE_GROUP)
# 0 - no; 1 - yes
# can set 15 samples for young/old and ventilator/not category



# no ventilator
samples <- subset(metadata, (SMTSD == "Lung") &
                    (DTHVNT == 0) &
                    (AGE_GROUP %in% c("20-29","30-39")))$SAMPID
tmp <- gtex[,sample(samples, size = min(15,length(samples)))]
saveRDS(tmp, "Outputs/3_GTExDataPrep/Ventilator/Lung_noventilator_young_sampled_data.rds")
rm(samples,tmp); gc()

samples <- subset(metadata, (SMTSD == "Lung") &
                    (DTHVNT == 0) &
                    (AGE_GROUP == "60-69"))$SAMPID
tmp <- gtex[,sample(samples, size = min(15,length(samples)))]
saveRDS(tmp, "Outputs/3_GTExDataPrep/Ventilator/Lung_noventilator_old_sampled_data.rds")
rm(samples,tmp); gc()

# ventilator
samples <- subset(metadata, (SMTSD == "Lung") &
                    (DTHVNT == 1) &
                    (AGE_GROUP %in% c("20-29","30-39")))$SAMPID
tmp <- gtex[,sample(samples, size = min(15,length(samples)))]
saveRDS(tmp, "Outputs/3_GTExDataPrep/Ventilator/Lung_ventilator_young_sampled_data.rds")
rm(samples,tmp); gc()

samples <- subset(metadata, (SMTSD == "Lung") &
                    (DTHVNT == 1) &
                    (AGE_GROUP == "60-69"))$SAMPID
tmp <- gtex[,sample(samples, size = min(15,length(samples)))]
saveRDS(tmp, "Outputs/3_GTExDataPrep/Ventilator/Lung_ventilator_old_sampled_data.rds")
rm(samples,tmp); gc()


# non-smoker
samples <- subset(metadata, (SMTSD == "Lung") &
                    (MHSMKSTS == "No") &
                    (AGE_GROUP %in% c("20-29","30-39")))$SAMPID
tmp <- gtex[,sample(samples, size = min(15,length(samples)))]
saveRDS(tmp, "Outputs/3_GTExDataPrep/Smoking/Lung_nonsmoker_young_sampled_data.rds")
rm(samples,tmp); gc()

samples <- subset(metadata, (SMTSD == "Lung") &
                    (MHSMKSTS == "No") &
                    (AGE_GROUP == "60-69"))$SAMPID
tmp <- gtex[,sample(samples, size = min(15,length(samples)))]
saveRDS(tmp, "Outputs/3_GTExDataPrep/Smoking/Lung_nonsmoker_old_sampled_data.rds")
rm(samples,tmp); gc()

# smoker
samples <- subset(metadata, (SMTSD == "Lung") &
                    (MHSMKSTS == "Yes") &
                    (AGE_GROUP %in% c("20-29","30-39")))$SAMPID
tmp <- gtex[,sample(samples, size = min(15,length(samples)))]
saveRDS(tmp, "Outputs/3_GTExDataPrep/Smoking/Lung_smoker_young_sampled_data.rds")
rm(samples,tmp); gc()

samples <- subset(metadata, (SMTSD == "Lung") &
                    (MHSMKSTS == "Yes") &
                    (AGE_GROUP == "60-69"))$SAMPID
tmp <- gtex[,sample(samples, size = min(15,length(samples)))]
saveRDS(tmp, "Outputs/3_GTExDataPrep/Smoking/Lung_smoker_old_sampled_data.rds")
rm(samples,tmp); gc()


# all smokers and all non-smokers, age-balanced
smk_samples <- c()
nonsmk_samples <- c()
for(age_group in c("20-29","30-39","40-49","50-59","60-69")){
  
  tmp <- subset(m, (MHSMKSTS == "No") & (AGE_GROUP == age_group))$SAMPID
  tmp <- sample(tmp, size = min(length(tmp),15), replace = F)
  nonsmk_samples <- c(nonsmk_samples, tmp)
  
  n <- length(tmp)
  rm(tmp)
  
  tmp <- subset(m, (MHSMKSTS == "Yes") & (AGE_GROUP == age_group))$SAMPID
  tmp <- sample(tmp, size = n, replace = F)
  smk_samples <- c(smk_samples, tmp)
  rm(tmp)
}

saveRDS(gtex[,nonsmk_samples], "Outputs/3_GTExDataPrep/Smoking/Lung_nonsmoker_sampled_data.rds")
saveRDS(gtex[,smk_samples], "Outputs/3_GTExDataPrep/Smoking/Lung_smoker_sampled_data.rds")

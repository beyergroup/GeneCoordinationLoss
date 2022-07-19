# Linear model of NMSE ~ Age + Expression

source("Scripts/functions.R")

net <- as.matrix(ReadRDS(paste0("Outputs/0_Preprocessing/stabsel_filtered_largestCC_network_Hs.rds")))

# retrieve ages from GTEx v7
metadata <- readRDS("/cellnet/AgeingClocks/Bulk/GTEx/Preprocessing/metadata.rds")
data <- readRDS("/cellnet/GeneCorrelation/Francisco/GTEx_analysis/04-Known_Covariates_Batch_Correction/Outputs/log2_batch_corrected_gtex.rds")
data <- t(data)
data <- data[,intersect(colnames(data),metadata$SAMPID)]
metadata <- metadata[match(colnames(data),metadata$SAMPID),]

tissues <- names(which(table(metadata$SMTSD) >= 350))
metadata <- subset(metadata, (SMTSD %in% tissues) | (SMTS == "Brain"))
metadata$Tissue <- as.character(metadata$SMTSD)
metadata$Tissue[metadata$SMTS == "Brain"] <- "Brain"
# metadata <- subset(metadata, AGE != "70-79")
data <- data[,as.character(metadata$SAMPID)]

n <- min(table(metadata$Tissue, metadata$AGE_GROUP)[,-6])

for(t in c(tissues,"Brain")){
  
  set.seed(seed = 1)
  ss <- c()
  for(a in c("20-29","30-39","40-49","50-59","60-69","70-79")){
    if(sum((metadata$Tissue == t) & (metadata$AGE_GROUP == a)) > n){
      ss <- c(ss, sample(as.character(subset(metadata, (Tissue == t) & (AGE_GROUP == a))$SAMPID),
                         size = n))
    } else{
      ss <- c(ss, as.character(subset(metadata, (Tissue == t) & (AGE_GROUP == a))$SAMPID))
    }
  }
  m <- subset(metadata, SAMPID %in% ss)
  d <- data[,as.character(m$SAMPID)]
  
  saveRDS(d, paste0("Outputs/5_Predictability/",gsub(" ","",t),"_v7_sampled_data.rds"))
  saveRDS(m, paste0("Outputs/5_Predictability/",gsub(" ","",t),"_v7_sampled_metadata.rds"))
  
  # center data
  centers <- rowMeans(d)
  d <- sweep(d, 1, centers, "-")
  
  saveRDS(d, paste0("Outputs/5_Predictability/",gsub(" ","",t),"_v7_sampled_centered_data.rds"))
  saveRDS(centers, paste0("Outputs/5_Predictability/",gsub(" ","",t),"_v7_sampled_centers.rds"))
  
  # network predictions
  d <- DataENSGToSymbol(d)
  pred <- PredictNet(net, as.matrix(d), maxiter = 1)
  saveRDS(pred, paste0("Outputs/5_Predictability/",gsub(" ","",t),"_v7_sampled_centered_net_predictions.rds"))
  
  # RMSE
  rmse <- sapply(intersect(rownames(pred),rownames(d)),
                 function(g) ComputeRMSEPerSample(pred = pred[g,], obs = d[g,]))
  rmse <- t(rmse)
  saveRDS(rmse, paste0("Outputs/5_Predictability/",gsub(" ","",t),"_v7_sampled_centered_rmse.rds"))
  
  # linear model RMSE ~ Age + Expression
  coefficients <- sapply(rownames(rmse), function(gene) {
    lm.data <- data.frame("RMSE" = rmse[gene,],
                          "Expression" = d[gene,],
                          "Age" = m$AGE)
    fit <- lm(RMSE ~ Expression + Age, data = lm.data)
    return(c("AgeSlope" = summary(fit)$coefficients["Age","Estimate"],
      "ExprSlope" = summary(fit)$coefficients["Expression","Estimate"],
      "AgePVal" = summary(fit)$coefficients["Age","Pr(>|t|)"],
      "ExprPVal" = summary(fit)$coefficients["Expression","Pr(>|t|)"]))
    # return(summary(fit)$coefficients["Age",c("Estimate","Pr(>|t|)")])
  })
  coefficients <- t(coefficients)
  saveRDS(coefficients, paste0("Outputs/5_Predictability/",gsub(" ","",t),
                               "_v7_sampled_coefficients.rds"))
  
  # # linear model NRMSE ~ Age
  # coefficients <- sapply(rownames(nrmse), function(gene) {
  #   lm.data <- data.frame("NRMSE" = nrmse[gene,],
  #                         "Age" = m$AGE)
  #   fit <- lm(NRMSE ~ Age, data = lm.data)
  #   return(c("AgeSlope" = summary(fit)$coefficients["Age","Estimate"],
  #            "AgePVal" = summary(fit)$coefficients["Age","Pr(>|t|)"]))
  # })
  # coefficients <- t(coefficients)
  # saveRDS(coefficients, paste0("Outputs/5_Predictability/",gsub(" ","",t),"_v7_sampled_coefficients_noexpr.rds"))
  
  # # regress out expression from NMSE
  # residuals <- sapply(rownames(nmse), function(gene) {
  #   lm.data <- data.frame("NMSE" = nmse[gene,],
  #                         "Expression" = d[gene,])
  #   fit <- lm(NMSE ~ Expression, data = lm.data)
  #   return(fit$residuals)
  # })
  # residuals <- t(residuals)
  # 
  # # feed residuals into limma
  # mm <- model.matrix(~ Age, data = data.frame("Age" = m$AGE))
  # fit <- limma::lmFit(residuals, mm)
  # fit <- eBayes(fit)
  
  rm(d,m,centers,pred,nmse,coefficients); gc()
}

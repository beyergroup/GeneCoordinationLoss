# Compute correlation to predictors

NET = "stabsel_filtered_trans_largestCC"
SLOPE_THRE = 0.005

source("Scripts/functions.R")

files <- list.files("Outputs/3_GTExDataPrep/Subset_Data/Max_Subset",
                    pattern = "centered_data", full.names = T)
files <- files[-grep("young|old",files)]

net <- ReadRDS(paste0("Outputs/0_Preprocessing/",NET,"_network_Hs.rds"))

tissues <- sapply(files, function(x) strsplit(tail(strsplit(x,"/")[[1]],1),
                                              "_")[[1]][1])

for(tissue in unique(tissues)){
  
  # read in files for different ages
  expr <- sapply(files[grep(tissue,files,fixed = T)],ReadRDS)
  expr <- lapply(expr, DataENSGToSymbol, remove_dup = T)
  names(expr) <- sapply(names(expr),
                        function(x) strsplit(tail(strsplit(x, split = "/")[[1]],1),
                                             "_")[[1]][2])
  
  mean_slopes <- data.frame("Gene" = NULL, "AvgSlope" = NULL)
  
  for(gene in intersect(rownames(expr[[1]]),rownames(net))){
    
    # reduce to predictive model
    prediction <- net[gene, net[gene,] != 0]
    
    if(!any(names(prediction) %in% rownames(expr[[1]])))
      next
    
    # compute correlation between gene and predictors
    cors <- lapply(expr,
                   function(m) apply(m[intersect(rownames(m),names(prediction)),,drop=F],
                                     1, cor, y = m[gene,]))
    cors <- do.call(cbind,cors)
    
    # multiply correlation coefficient by sign of the relationship between genes
    cors <- cors*sign(prediction)[rownames(cors)]
    
    # slope of correlations with age (70s)
    lm_results <- t(apply(cors, 1, function(cor){
      fit <- lm(data = data.frame("Correlation" = cor,
                                  "Age" = c(25,35,45,55,65,70)[1:length(cor)]),
                formula = Correlation ~ Age)
      coefs <- summary(fit)$coefficients
      return(c("Slope" = coefs["Age","Estimate"],
               "pval" = coefs["Age","Pr(>|t|)"]))
    }))
      
    # slope of correlations with age (no 70s)
    lm_results_no70 <- t(apply(cors, 1, function(cor){
      fit <- lm(data = data.frame("Correlation" = cor[1:5],
                                  "Age" = c(25,35,45,55,65)),
                formula = Correlation ~ Age)
      coefs <- summary(fit)$coefficients
      return(c("Slope" = coefs["Age","Estimate"],
               "pval" = coefs["Age","Pr(>|t|)"]))
    }))
    
    # sig <- names(which((lm_results[,"pval"] < 0.1) &
    #                      (lm_results_no70[,"pval"] < 0.1)))
    
    # look at slope effect size rather than p-values
    sig <- names(which(((lm_results[,"Slope"] < -SLOPE_THRE) &
                          (lm_results_no70[,"Slope"] < -SLOPE_THRE)) |
                         ((lm_results[,"Slope"] > SLOPE_THRE) &
                            (lm_results_no70[,"Slope"] > SLOPE_THRE))))
    
    if(length(sig) > 0){
      
      # # look for consistent slope signs
      # sig <- sig[sign(lm_results[sig,"Slope"]) ==
      #              sign(lm_results_no70[sig,"Slope"])]
      
      # get general trend of slope by weighted average with module of relationship weights
      mean_slopes <- rbind.data.frame(mean_slopes,
                           data.frame("Gene" = gene,
                                      "AvgSlope" = mean(lm_results[sig,"Slope"],
                                                        weights = abs(net[gene,sig]))))
    }
  }
  
  WriteRDS(mean_slopes, paste0("Outputs/5_Predictability/Age/Age_MaxSubset/",
                               tissue,"_predcor_robust_slope.rds"))
}

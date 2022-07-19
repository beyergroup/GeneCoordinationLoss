# Correlation between pred-obs cor. and age

args = commandArgs(trailingOnly=TRUE)
NET = args[1]


files <- list.files(paste0("Outputs/Human_Network/",NET,"/Predictability/AgeTissue"),
                    pattern = "sampled_correlations.rds", full.names = T)
tissues <- unique(sapply(files, function(x) strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1]))

ages <- c(25,35,45,55,65)


for(tissue in tissues){
  
  cors <- sapply(files[grep(tissue,files,fixed = T)], readRDS)
  age_cors <- t(apply(cors, 1,
                      function(x) tryCatch(unlist(cor.test(x, y = ages,
                                                           method = "pearson")[c("estimate","p.value")]),
                                           error = function(e) return(c(NA,NA)))))
  colnames(age_cors) <- c("cor","pval")
  age_cors <- cbind(age_cors, "FDR" = p.adjust(age_cors[,"pval"]))
  
  saveRDS(age_cors,
          paste0("Outputs/Human_Network/",NET,"/Predictability/AgeTissue/",tissue,"_sampled_age_corP.rds"))
}

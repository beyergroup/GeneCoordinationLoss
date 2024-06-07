
SLOPE_DIR = "Outputs/5_Predictability/Age/Age_MaxSubset"
COR_THRE = "well_predicted"

source("Scripts/functions.R")


files <- list.files(SLOPE_DIR, full.names = T,
                    pattern = paste0("ageslope_",COR_THRE,"_pearson"))

pdf(paste0("Plots/5_Predictability/Age/Age_MaxSubset/compare_pearson_spearman_",
           COR_THRE,".pdf"))

for(file in files){
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  pearson <- ReadRDS(file)[,"Slope"]
  spearman <- ReadRDS(gsub("_pearson.rds",".rds",file))[,"Slope"]
  
  plot.data <- data.frame("Spearman" = spearman[intersect(names(pearson),names(spearman))],
                          "Pearson" = pearson[intersect(names(pearson),names(spearman))])
  
  print(ggplot(plot.data) +
    geom_point(aes(x = Pearson, y = Spearman)) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Pearson-based slopes") + ylab("Spearman-based slopes") +
    ggtitle(tissue, paste("Correlation threshold:",COR_THRE)) +
    theme_classic() + theme(text = element_text(size = 20)))
}

dev.off()

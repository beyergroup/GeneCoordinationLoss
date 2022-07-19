.libPaths("Resources/Rlibs/R-4.0.3")
source("Scripts/functions.R")
library(ggplot2)

files <- list.files("Outputs/5_Predictability",
                    pattern = "correlations")
# files <- files[grep("young|old",files)]
files <- files[-grep("young|old",files)]

df <- data.frame()

for(mode in c("cis","trans")){
  
  f <- files[grep(mode,files)]
  
  for(file in f){
    
    tissue <- strsplit(file,"_")[[1]][1]
    cor <- readRDS(paste0("Outputs/5_Predictability/",file))
    
    df <- rbind.data.frame(df,
                           data.frame("Correlation" = cor,
                                      "Tissue" = tissue,
                                      "Mode" = mode))
  }
}

for(file in gsub("_cis","",files[grep("cis",files)])){
  
  tissue <- strsplit(file,"_")[[1]][1]
  cor <- readRDS(paste0("Outputs/5_Predictability/WellPredicted_TissueFilters/",file))
  
  df <- rbind.data.frame(df,
                         data.frame("Correlation" = cor,
                                    "Tissue" = tissue,
                                    "Mode" = "global"))
}

for(file in gsub("_cis","",files[grep("cis",files)])){
  
  tissue <- strsplit(file,"_")[[1]][1]
  cor <- readRDS(paste0("Outputs/5_Predictability/WellPredicted_TissueFilters/Randomized/",file))
  
  df <- rbind.data.frame(df,
                         data.frame("Correlation" = cor,
                                    "Tissue" = tissue,
                                    "Mode" = "random"))
}


pdf("Plots/5_Predictability/cis_trans_tissues.pdf",
    width = 16, height = 8)
ggplot(df) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_density(aes(x = Correlation, color = Mode)) +
  scale_color_manual(values = c("cis" = cis_color, "trans" = trans_color,
                                "global" = net_color, "random" = "grey"),
                     labels = c("cis" = "Cis-predictors only",
                                "trans" = "Trans-predictors only",
                                "global" = "Full",
                                "random" = "Randomized")) +
  labs(color = "Model") + xlab("Correlation predicted - observed") +
  facet_wrap(~ Tissue) + theme(text = element_text(size = 20),
                               legend.position = "bottom",
                               axis.title.y = element_blank())
dev.off()

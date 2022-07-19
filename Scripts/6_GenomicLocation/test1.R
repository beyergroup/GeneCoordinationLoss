
.libPaths("Resources/Rlibs/R-4.0.3")
source("Scripts/functions.R")
library(ggplot2)

files <- list.files("Outputs/5_Predictability",
                    pattern = "correlations_spearman")
files <- files[grep("young|old",files)]

df <- data.frame()

for(mode in c("cis","trans")){
  
  f <- files[grep(mode,files)]
  
  for(file in f){
    
    tissue <- strsplit(file,"_")[[1]][2]
    cor <- readRDS(paste0("Outputs/5_Predictability/",file))
    
    df <- rbind.data.frame(df,
                           data.frame("Correlation" = cor,
                                      "Tissue" = tissue, 
                                      "Age" = strsplit(file,"_")[[1]][1],
                                      "Mode" = mode))
  }
}

for(file in files[-c(grep("cis",files),grep("trans",files))]){
  
  tissue <- strsplit(file,"_")[[1]][2]
  cor <- readRDS(paste0("Outputs/5_Predictability/",file))
  
  df <- rbind.data.frame(df,
                         data.frame("Correlation" = cor,
                                    "Tissue" = tissue,
                                    "Age" = strsplit(file,"_")[[1]][1],
                                    "Mode" = "global"))
}

pdf("Plots/5_Predictability/cis_trans_young_old_spearman.pdf",
    width = 16, height = 8)
ggplot(df) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_density(aes(x = Correlation, color = Mode, linetype = Age)) +
  scale_color_manual(values = c("cis" = cis_color, "trans" = trans_color,
                                "global" = net_color),
                     labels = c("cis" = "Cis-predictors only",
                                "trans" = "Trans-predictors only",
                                "global" = "Full model")) +
  scale_linetype_manual(values = c("young" = "solid", "old" = "dashed"),
                        labels = c("young" = "Young", "old" = "Old")) +
  labs(color = "Model") + xlab("Correlation (Spearman) predicted - observed") +
  facet_wrap(~ Tissue) + theme(text = element_text(size = 20),
                               legend.position = "bottom",
                               axis.title.y = element_blank())
dev.off()

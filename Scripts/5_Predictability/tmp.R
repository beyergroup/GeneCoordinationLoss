files <- list.files("Outputs/5_Predictability/tmp",
                    pattern = paste0("ageslope_",COR_THRE),
                    full.names = T)

plot.data <- data.frame()

for(file in files){
  
  slopes <- ReadRDS(file)
  background <- ReadRDS(gsub("ageslope","ageslope_background",file))
  
  plot.data <- rbind.data.frame(plot.data,
                                data.frame("Slope" = slopes[,"AdjustedSlope"],
                                           "BG" = background[,"AdjustedSlope"],
                                           "Tissue" = strsplit(tail(strsplit(file,"/")[[1]],
                                                                    1),
                                                               "_")[[1]][1]))
  
  rm(slopes,backgrounds); gc()
}

pdf(paste0("Plots/5_Predictability/tmp/slopes_LFCs_",COR_THRE,".pdf"),
    width = 11, height = 11)
ggplot(plot.data) +
  geom_density(aes(x = BG), color = "dimgrey") +
  geom_density(aes(x = Slope), color = "red", linetype = "dashed") +
  facet_wrap(~ Tissue, ncol = 3) +
  xlab("Slope distribution") +
  ggtitle(paste0("Min avg cor = ",COR_THRE)) +
  theme(text = element_text(size = 20),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

rm(plot.data); gc()


# Call hits based on tails of background dists

hit_list <- list()

for(file in files){
  
  tissue <- strsplit(tail(strsplit(file,"/")[[1]],1),"_")[[1]][1]
  
  background <- ReadRDS(gsub("ageslope","ageslope_background",file))
  
  thre <- quantile(background[,"AdjustedSlope"], probs = c(0.025, 0.975))
  
  slopes <- ReadRDS(file)
  pos <- rownames(slopes)[slopes[,"AdjustedSlope"] > thre[2]]
  neg <- rownames(slopes)[slopes[,"AdjustedSlope"] < thre[1]]
  hit_list[[tissue]] <- list("PositiveSlope" = pos,
                             "NegativeSlope" = neg,
                             "Unchanged" = setdiff(rownames(slopes),c(pos,neg)))
  
  rm(slopes,background,tissue,pos,neg); gc()
}

WriteRDS(hit_list, paste0("Outputs/5_Predictability/tmp/hit_list_slopes_",
                          COR_THRE,".rds"))

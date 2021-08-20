# Scatterplot of mean predictability across tissues and mean absolute
# expression levels for each module

library(ggplot2, lib.loc = "Resources/Rlibs/R-4.0.3/")

for(NET in c("stabsel","stabsel_pcclasso")){
  
  for(METHOD in c("greedy","rwalk")){
    
    m <- readRDS(list.files(paste0("Outputs/Human_Network/",NET,"/Topology/Modules"),
      full.names = T, recursive = T, pattern = paste("mean_predictability", METHOD, sep = "_")))
    
    abs_expr <- readRDS(paste0("Outputs/Human_Network/",NET,
                               "/Topology/Modules/TissueDE/CrossTissue_mean_absexpression_",
                               METHOD,".rds"))
    abs_expr <- abs_expr[rownames(m)]
    
    plot.data <- data.frame("CrossTissuePredictability" = m[,"CrossTissue"],
                            "MeanTissuePredictability" = rowMeans(m[,colnames(m) != "CrossTissue"]),
                            "AbsoluteExpression" = abs_expr)
    pdf(paste0("Plots/Human_Network/",NET,"/Predictability/modules_",
               METHOD,"_predictabilityvsexpr.pdf"))
    print(ggplot(plot.data) +
      geom_point(aes(x = AbsoluteExpression, y = CrossTissuePredictability)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      xlab("Normalized cross-tissue expression")+
      ylab("Cross-tissue predictability") +
      theme_classic() + theme(text = element_text(size = 20)))
    print(ggplot(plot.data) +
      geom_point(aes(x = AbsoluteExpression, y = MeanTissuePredictability)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      xlab("Normalized cross-tissue expression")+
      ylab("Mean predictability") +
      theme_classic() + theme(text = element_text(size = 20)))
    dev.off()
  }
  
}

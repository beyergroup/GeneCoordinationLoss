# Plot predictability

TYPE = "Tabula_Sapiens"
exclude_poorly_predicted = T
net = "stabsel"
method = "smartseq2"

correlations <- list()

files <- list.files(paste0("Outputs/Human_Network/",net,"/Predictability/",TYPE),
                    pattern = "correlations.rds", full.names = T)
files <- files[grep(method,files)]

for(file in files){
  
  prefix <- strsplit(tail(strsplit(file, "/")[[1]],1),"_sampled_")[[1]][1]
  
  correlations[[prefix]] <- readRDS(file)
  
}


genes <- lapply(correlations, names)
intersect <- c()
for(gene in genes[[which.min(lapply(genes,length))]]){
  if(all(unlist(lapply(genes, function(gene_vector) gene %in% gene_vector))))
    intersect <- c(intersect,gene)
}
if(exclude_poorly_predicted){
  poorly_predicted <- readRDS(paste0("Outputs/Human_Network/",net,
                                     "/Predictability/Tissue/poorly_predicted_crosstissue.rds"))
  intersect <- intersect[!(intersect %in% poorly_predicted)]
}
correlations <- lapply(correlations, function(x) x[intersect])

length(intersect)


library(reshape2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggplot2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggpubr, lib.loc = "Resources/Rlibs/R-4.0.3/")

plot.data <- reshape2::melt(correlations)
plot.data$Tissue <- sapply(plot.data$L1, function(x) strsplit(x,"_")[[1]][1])
plot.data$Method <- sapply(plot.data$L1, function(x) tail(strsplit(x,"_")[[1]],1))
plot.data$CellType <- sapply(1:nrow(plot.data),
                             function(i) paste(strsplit(plot.data$L1[i],
                                                        "_")[[1]][2:(grep(plot.data$Method[i],
                                                                         strsplit(plot.data$L1[i],"_")[[1]])-1)],
                                               collapse = " "))
plot.data$Method[plot.data$Method == "smartseq2"] <- "SS2"

plots <- list()
for(d in unique(plot.data$L1)){
  tissue <- unique(subset(plot.data, L1 == d)$Tissue)
  celltype <- unique(subset(plot.data, L1 == d)$CellType)
  plots[[d]] <- ggplot(subset(plot.data, L1 == d)) +
    geom_density(aes(x = value), fill = "blue", alpha = .4, size = 1) +
    xlim(c(-1,1)) + geom_vline(xintercept = 0, linetype = "dashed") +
    xlab("Correlation coefficient") +
    ggtitle(celltype, subtitle = paste0(tissue," (",
                                        c("10X" = "10X", "smartseq2" = "SS2")[method],")")) +
    theme_classic() + theme(text = element_text(size = 20),
                            legend.position = "bottom", legend.title = element_blank(),
                            axis.title.y = element_blank(), title = element_text(size = 18))
}

dims <- list("smartseq2" = c("height" = 5, "width" = 10),
             "10X" = c("height" = 9, "width" = 13))

pdf(paste0("Plots/Human_Network/stabsel/Predictability/predictability_Tabula_Sapiens_tmp_",method,".pdf"),
    height = dims[[method]][1], width = dims[[method]][2])
ggarrange(plotlist = plots, align = "hv", common.legend = T, legend = "bottom")
dev.off()

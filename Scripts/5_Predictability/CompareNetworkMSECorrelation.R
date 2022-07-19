.libPaths("Resources/Rlibs/R-4.0.3/")
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
source("Scripts/functions.R")


mse_files <- list.files("Outputs/5_Predictability", pattern = "_sampled_centered_NMSE.rds", full.names = T)
cor_files <- list.files("Outputs/5_Predictability", pattern = "_sampled_centered_correlations.rds", full.names = T)


well_predicted <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")


# Scatterplot

pdf("Plots/5_Predictability/correlation_NMSE.pdf")

for(i in 1:length(mse_files)){
  
  mse <- readRDS(mse_files[i])
  cor <- readRDS(cor_files[i])
  
  tissue <- strsplit(tail(strsplit(mse_files[i], split = "/")[[1]],1), "_")[[1]][1]
  age_group <- strsplit(tail(strsplit(mse_files[i], split = "/")[[1]],1), "_")[[1]][2]
  
  genes <- intersect(names(cor),names(mse))
  genes <- intersect(genes, well_predicted[[tissue]])
  
  print(smoothScatter(x = log10(mse[genes]), y = cor[genes],
                      main = paste0(tissue,"\n",age_group),
                      xlab = "log10 NMSE", ylab = "Correlation",
                      ylim = c(-1,1)))
}

dev.off()



cor_files <- list.files("Outputs/5_Predictability", pattern = "corLFCs_YvsO.rds", full.names = T)
delta_files <- list.files(paste0("Outputs/Human_Network/stabsel/Predictability/AgeTissue"),
                        pattern = "_sampled_ageDP.rds", full.names = T)

tissues <- sapply(cor_files, function(x) strsplit(tail(strsplit(x,"/")[[1]],1),"_")[[1]][1])

plot.data <- data.frame()

for(tissue in tissues){
  
  corLFCs <- readRDS(cor_files[grep(tissue,cor_files,fixed=T)])
  deltaLFCs <- readRDS(delta_files[grep(tissue,delta_files,fixed=T)])
  deltaLFCs <- limma::topTable(fit = deltaLFCs, coef = "AgeGroupOld", number = nrow(deltaLFCs))
  
  common <- intersect(rownames(deltaLFCs),names(corLFCs))
  common <- intersect(common,well_predicted[[tissue]])
  
  plot.data <- rbind.data.frame(plot.data,
                                data.frame("DeltaLFC" = deltaLFCs[common,"logFC"],
                                           "CorLFC" = corLFCs[common],
                                           "Tissue" = tissue,
                                           "Gene" = common))
  rm(corLFCs,deltaLFCs,common,tissue); gc()
}

ggplot(plot.data) +
  geom_point(aes(x = DeltaLFC, y = CorLFC)) +
  facet_wrap(~ Tissue) +
  coord_cartesian(ylim = c(-5,5))

# overlap is total shit (they're actually anti-correlated wtf)



# Look at top cases (Lung, Skin and Artery)

gene = "MICAL2"
gene = "ADAM30"

pdf(paste0("Plots/5_Predictability/",gene,"_corLFCs.pdf"),
    width = 5, height = 5)
# for(tissue in tissues){
for(tissue in c("Esophagus-Mucosa","Lung")){
# for(tissue in c("Skin-NotSunExposed(Suprapubic)","Thyroid")){
  
  tmp <- subset(plot.data, Tissue == tissue)
  tmp <- tmp[order(abs(tmp$CorLFC), decreasing = T),]
  genes <- head(tmp$Gene,5)
  
  observed_young <- readRDS(paste0("Outputs/5_Predictability/young_",tissue,"_sampled_centered_data.rds"))
  observed_old <- readRDS(paste0("Outputs/5_Predictability/old_",tissue,"_sampled_centered_data.rds"))
  observed_young <- DataENSGToSymbol(observed_young)
  observed_old <- DataENSGToSymbol(observed_old)
  predicted_young <- readRDS(paste0("Outputs/5_Predictability/young_",tissue,"_sampled_centered_net_predictions.rds"))
  predicted_old <- readRDS(paste0("Outputs/5_Predictability/old_",tissue,"_sampled_centered_net_predictions.rds"))
  
  # centers_young <- readRDS(paste0("Outputs/5_Predictability/young_",tissue,"_sampled_centers.rds"))
  # names(centers_young) <- VectorENSGToSymbol(names(centers_young))
  # centers_old <- readRDS(paste0("Outputs/5_Predictability/old_",tissue,"_sampled_centers.rds"))
  # names(centers_old) <- VectorENSGToSymbol(names(centers_old))
  
  # pdf(paste0("Plots/5_Predictability/",tissue,"_top_corLFCs.pdf"),
  #     width = 5, height = 5)
  # for(gene in genes){
    
    p.data <- rbind.data.frame(data.frame("Observed" = observed_young[gene,colnames(observed_young)],
                                          "Predicted" = predicted_young[gene,colnames(observed_young)],
                                          "Age" = "Young"),
                                 data.frame("Observed" = observed_old[gene,colnames(observed_old)],
                                          "Predicted" = predicted_old[gene,colnames(predicted_old)],
                                          "Age" = "Old"))
    p.data <- melt(p.data, measure.vars = c("Observed","Predicted"),
                   variable.name = "Type")
    p.data$Age <- factor(as.character(p.data$Age), levels = c("Young","Old"))
    
    if(grepl("-",tissue)){
      tissue_title <- gsub("-"," - ",tissue)
    } else{
      tissue_title <- tissue
    }
    
    print(ggplot(p.data) +
            geom_violin(aes(x = Age, y = value, fill = Type),
                        draw_quantiles = 0.5, size = 1.5) +
            geom_point(aes(x = Age, y = value, color = Type),
                       position = position_jitterdodge(dodge.width = .8)) +
            scale_color_manual(values = c("Observed" = "black", "Predicted" = "black")) +
            scale_fill_manual(values = c("Observed" = obs_color, "Predicted" = net_color)) +
            ylab(paste0(gene, " centered expr.")) +
            # ylab(paste0(gene, " expression")) +
            ggtitle(tissue_title, subtitle = paste0("Error LFC = ",
                                              round(tmp[gene,"CorLFC"],2))) +
            theme_classic() + theme(text = element_text(size = 20),
                                    legend.title = element_blank(),
                                    legend.position = "bottom"))
  # }
  # dev.off()
  
}
dev.off()

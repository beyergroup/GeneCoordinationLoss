# For each tissue, heatmap of absolute expression, observed FC and predicted FC
# of genes with significant age-related changes in predictability

source("Scripts/functions.R")
library(pheatmap, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggpubr, lib.loc = "Resources/Rlibs/R-4.0.3/")

for(NET in c("stabsel","stabsel_pcclasso")){
  
  dir.create(paste0("Plots/Human_Network/",NET,"/Predictability/Age"))
  
  files <- list.files(paste0("Outputs/Human_Network/",NET,
                             "/Predictability/AgeTissue"),
                      pattern = "sampled_ageDP.rds", full.names = T)
  
  for(file in files){
    
    ageDP <- readRDS(file)
    ageDP <- limma::topTable(fit = ageDP, coef = "AgeGroupOld", number = nrow(ageDP))
    
    if(sum(ageDP$adj.P.Val < 0.05) > 1){
      
      tissue <- head(strsplit(tail(strsplit(file, "/")[[1]],1), "_")[[1]],1)
      ageDP <- subset(ageDP, adj.P.Val < 0.05)
      genes <- rownames(ageDP)[order(ageDP$logFC, decreasing = T)]
      
      Age_ann <- c(rep("20-29",31),rep("30-39",31),rep("50-59",31),rep("60-69",31))
      
      # centered expression
      centered <- cbind(readRDS(paste0("GTEx_Networks/AgeTissue_Networks/Outputs/",tissue,
                                       "_20-29_sampled_centered_data.rds")),
                        readRDS(paste0("GTEx_Networks/AgeTissue_Networks/Outputs/",tissue,
                                       "_30-39_sampled_centered_data.rds")),
                        readRDS(paste0("GTEx_Networks/AgeTissue_Networks/Outputs/",tissue,
                                       "_50-59_sampled_centered_data.rds")),
                        readRDS(paste0("GTEx_Networks/AgeTissue_Networks/Outputs/",tissue,
                                       "_60-69_sampled_centered_data.rds")))
      centered <- DataENSGToSymbol(centered, remove_dup = T)
      centered <- centered[genes, , drop = F]
      names(Age_ann) <- colnames(centered)
      
      pheatmap(centered, scale = "none", cluster_rows = F, cluster_cols = F,
               breaks = seq(from = -(max(abs(preds))), to = max(abs(preds)), length.out = 501),
               color = colorRampPalette(c("blue","white","red"))(n = 501),
               annotation_col = data.frame("Age" = as.factor(Age_ann)),
               annotation_colors = list("Age" = age_palette),
               show_colnames = F, gaps_col = cumsum(rep(31,4)),
               cellwidth = 2, cellheight = 10,
               main = paste0(tissue," centered expression"),
               filename = paste0("Plots/Human_Network/",NET,"/Predictability/Age/",
                                 tissue,"_centered_heatmap.pdf"))
      
      
      # predictions
      preds <- cbind(readRDS(paste0("Outputs/Human_Network/",NET,"/Predictability/AgeTissue/",tissue,
                                    "_20-29_sampled_net_predictions.rds")),
                     readRDS(paste0("Outputs/Human_Network/",NET,"/Predictability/AgeTissue/",tissue,
                                    "_30-39_sampled_net_predictions.rds")),
                     readRDS(paste0("Outputs/Human_Network/",NET,"/Predictability/AgeTissue/",tissue,
                                    "_50-59_sampled_net_predictions.rds")),
                     readRDS(paste0("Outputs/Human_Network/",NET,"/Predictability/AgeTissue/",tissue,
                                    "_60-69_sampled_net_predictions.rds")))
      preds <- preds[genes, , drop = F]
      
      pheatmap(preds, scale = "none", cluster_rows = F, cluster_cols = F,
               breaks = seq(from = -(max(abs(preds))), to = max(abs(preds)), length.out = 501),
               color = colorRampPalette(c("blue","white","red"))(n = 501),
               annotation_col = data.frame("Age" = as.factor(Age_ann)),
               annotation_colors = list("Age" = age_palette),
               show_colnames = F, gaps_col = cumsum(rep(31,4)),
               cellwidth = 2, cellheight = 10,
               main = paste0(tissue," predictions"),
               filename = paste0("Plots/Human_Network/",NET,"/Predictability/Age/",
                                 tissue,"_predictions_heatmap.pdf"))
      
      
      # absolute expression
      abs_expr <- readRDS("Data/GTEx/DESeq2_normalized_gtex.rds")
      abs_expr <- abs_expr[,colnames(preds)]
      abs_expr <- DataENSGToSymbol(abs_expr, remove_dup = T)
      abs_expr <- log2(1+abs_expr[genes, , drop = F])
      
      pheatmap(abs_expr, scale = "none", cluster_rows = F, cluster_cols = F,
               color = colorRampPalette(c("white","red"))(n = 501),
               annotation_col = data.frame("Age" = as.factor(Age_ann)),
               annotation_colors = list("Age" = age_palette),
               show_colnames = F, gaps_col = cumsum(rep(31,4)),
               cellwidth = 2, cellheight = 10,
               main = paste0(tissue," norm. expression"),
               filename = paste0("Plots/Human_Network/",NET,"/Predictability/Age/",
                                 tissue,"_absexpression_heatmap.pdf"))
      
      
      # deviation
      delta <- cbind(readRDS(paste0("Outputs/Human_Network/",NET,"/Predictability/AgeTissue/",tissue,
                                    "_20-29_sampled_delta.rds")),
                     readRDS(paste0("Outputs/Human_Network/",NET,"/Predictability/AgeTissue/",tissue,
                                    "_30-39_sampled_delta.rds")),
                     readRDS(paste0("Outputs/Human_Network/",NET,"/Predictability/AgeTissue/",tissue,
                                    "_50-59_sampled_delta.rds")),
                     readRDS(paste0("Outputs/Human_Network/",NET,"/Predictability/AgeTissue/",tissue,
                                    "_60-69_sampled_delta.rds")))
      delta <- delta[genes, , drop = F]
      
      pheatmap(delta, scale = "none", cluster_rows = F, cluster_cols = F,
               color = colorRampPalette(c("white","red"))(n = 501),
               annotation_col = data.frame("Age" = as.factor(Age_ann)),
               annotation_colors = list("Age" = age_palette),
               show_colnames = F, gaps_col = cumsum(rep(31,4)),
               cellwidth = 2, cellheight = 10,
               main = paste0(tissue," deviations"),
               filename = paste0("Plots/Human_Network/",NET,"/Predictability/Age/",
                                 tissue,"_delta_heatmap.pdf"))
    }
  }
}

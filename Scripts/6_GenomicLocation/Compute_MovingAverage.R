# Moving average of LFCs

source("Scripts/functions.R")


# Well-predicted genes
well_predicted_genes <- readRDS("Outputs/5_Predictability/WellPredicted_TissueFilters/well_predicted_genes.rds")

# Gene data
all.genes <- ReadRDS("Outputs/6_GenomicLocation/gene_info.rds")


moving_averages <- list()

for(tissue in names(well_predicted_genes)[-4]){
  
  ageLFCs <- ReadRDS(paste0("Outputs/5_Predictability/",tissue,"_trans_corLFCs_YvsO_adj_spearman.rds"))
  targets <- intersect(names(ageLFCs),all.genes$hgnc_symbol)
  
  data <- data.frame("LFC" = ageLFCs[targets],
                     "Gene" = targets,
                     "Position" = all.genes[targets,"Position"],
                     "Chromosome" = all.genes[targets,"chromosome_name"],
                     "Tissue" = tissue)
  
  for(WINDOW in c(1,2,3)){
    
    data[,paste0("MA",WINDOW)] <- NA
    
    data <- data[order(data$Position),]
    for(i in (WINDOW+1):(nrow(data)-WINDOW)){
      
      # check that window is not crossing chromosome boundaries
      if(data[(i-WINDOW),"Chromosome"] == data[i,"Chromosome"] &
         (data[(WINDOW+i),"Chromosome"] == data[i,"Chromosome"])){
        
        data[i,paste0("MA",WINDOW)] <- mean(data[(i-WINDOW):(i+WINDOW),"LFC"])
      }
    }
  }
  
  moving_averages[[tissue]] <- data
  rm(ageLFCs, targets, data); gc()
}

WriteRDS(moving_averages, "Outputs/6_GenomicLocation/moving_averages_trans.rds")

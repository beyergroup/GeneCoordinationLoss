# Count sample numbers in given subsetting technique for Supplemantary Tables

args <- commandArgs(trailingOnly = TRUE)
WDIR = args[1]
OUTPUT = args[2]

source("Scripts/functions.R")
library(reshape)
library(openxlsx)
source("Scripts/3_GTExDataPrep/params.R")

files <- list.files(WDIR, pattern = "sampled_data.rds")

df <- data.frame("Tissue" = NULL, "Age" = NULL, "Samples" = NULL)

for(file in files){
  
  d <- ReadRDS(paste0(WDIR,"/",file))
  t <- strsplit(file,"_")[[1]][1]
  a <- strsplit(file,"_")[[1]][2]
  
  if(!(t %in% RESTRICTED_TISSUES))
    next
  
  df <- rbind.data.frame(df,
                         data.frame(Tissue = t, Age = a, Samples = ncol(d)))
  
  rm(d,t,a)
}

# Place in correct format
table <- cast(df, Tissue ~ Age)

# Write to excel file
write.xlsx(table, file = OUTPUT)

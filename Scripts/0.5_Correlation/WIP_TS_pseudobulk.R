library(reshape2)
library(ggplot2)
source("Scripts/functions.R")

core.complex.list <- ReadRDS("Outputs/0.5_Correlation/coreComplexes.rds")

files_SS2 <- list.files("Outputs/3_TSDataPrep/smartseq2/Pseudobulk/Subsets")
files_SS2 <- files_SS2[-grep("50.rds",files_SS2)]
files_10X <- list.files("Outputs/3_TSDataPrep/10X/Pseudobulk/Subsets")

# ------------------------------------------------------------------------------

# Get DESeq2-normalized Tabula Sapiens for expression and variance computation

abs_expression <- list("Smart-seq2" = list(), "10X" = list())
variance <- list("Smart-seq2" = list(), "10X" = list())
within_cor <- list("Smart-seq2" = list(), "10X" = list())

for(file in files_SS2){
  
  cell_type <- gsub(".rds","",tail(strsplit(file,"_")[[1]],1))
  
  data <- ReadRDS(paste0("Outputs/3_TSDataPrep/smartseq2/Pseudobulk/Subsets/",
                         file))
  data <- log2(data+1)
  
  abs_expression[["Smart-seq2"]][[cell_type]] <- rowMeans(data)
  variance[["Smart-seq2"]][[cell_type]] <- apply(data, 1, var, na.rm = T)
  
  within_cor[["Smart-seq2"]][[cell_type]] <- lapply(core.complex.list,
                                                    function(genes){
                                                      d <- data[intersect(genes,rownames(data)),]
                                                      d[d == 0] <- NA
                                                      return(cor(t(d), use = "pairwise.complete.obs"))})
}

for(file in files_10X){
  
  cell_type <- gsub(".rds","",tail(strsplit(file,"_")[[1]],1))
  
  data <- ReadRDS(paste0("Outputs/3_TSDataPrep/10X/Pseudobulk/Subsets/",
                         file))
  data <- log2(data+1)
  
  abs_expression[["10X"]][[cell_type]] <- rowMeans(data)
  variance[["10X"]][[cell_type]] <- apply(data, 1, var, na.rm = T)
  
  within_cor[["10X"]][[cell_type]] <- lapply(core.complex.list,
                                             function(genes){
                                               d <- data[intersect(genes,rownames(data)),]
                                               # d <- d[rowSums(d == 0) == 0,]
                                               # if(length(d) > 10){
                                               #   return(cor(t(d)))
                                               # } else{
                                               #   return(NA)
                                               # }
                                               d[d == 0] <- NA
                                               return(cor(t(d), use = "pairwise.complete.obs"))})
}

saveRDS(abs_expression, "Outputs/0.5_Correlation/avg_absexpression_TS.rds")
saveRDS(variance, "Outputs/0.5_Correlation/expressionvar_TS.rds")
saveRDS(within_cor, "Outputs/0.5_Correlation/withincor_corecomplexes_TS.rds")

# ------------------------------------------------------------------------------

# Overall correlation dist of core complexes

plot.data <- lapply(within_cor,
                    function(l) lapply(l,
                                       function(ll) lapply(ll, function(mat) {
                                         if(ncol(mat) > 2){
                                           return(mean(mat[upper.tri(mat)],
                                                       na.rm = TRUE))
                                         } else{
                                           return(NA)}})))

plot.data <- rbind.data.frame(cbind.data.frame(melt(plot.data$`Smart-seq2`,
                                                    value.name = "Correlation"),
                                               "Method" = "Smart-seq2"),
                              cbind.data.frame(melt(plot.data$`10X`,
                                                    value.name = "Correlation"),
                                               "Method" = "10X"))
plot.data <- plot.data[!is.na(plot.data$Correlation),]
plot.data <- droplevels.data.frame(plot.data)

pdf("Plots/0.5_Correlation/corecomplex_meancor_TS.pdf", width = 55, height = 25)
ggplot(plot.data) +
  geom_point(aes(x = L2, y = Correlation, color = L1, shape = Method)) +
  scale_color_manual(values = cell_type_palette, name = "Cell Type") +
  xlab("Core complexes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 20), legend.position = "top")
ggplot(subset(plot.data, Method == "Smart-seq2")) +
  geom_point(aes(x = L2, y = Correlation, color = L1)) +
  scale_color_manual(values = cell_type_palette, name = "Cell Type") +
  xlab("Core complexes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 20), legend.position = "top")
dev.off()



tmp <- subset(plot.data, Method == "Smart-seq2")
whee <- c()
for(complex in unique(tmp$L2)){
  whee <- c(whee, min(subset(tmp, L2 == complex)$Correlation))
}
names(whee) <- unique(tmp$L2)

# ------------------------------------------------------------------------------

# Plot average complex expression

plot.data <- data.frame()

for(c in names(core.complex.list)){
  
  avg <- rbind.data.frame(cbind.data.frame(melt(lapply(abs_expression[["10X"]],
                                                       function(x) mean(x[intersect(names(x),
                                                                                    core.complex.list[[c]])])),
                                                value.name = "Mean"),
                                           Method = "10X"),
                          cbind.data.frame(melt(lapply(abs_expression[["Smart-seq2"]],
                                                       function(x) mean(x[intersect(names(x),
                                                                                    core.complex.list[[c]])])),
                                                value.name = "Mean"),
                                           Method = "Smart-seq2"))
  var <- rbind.data.frame(cbind.data.frame(melt(lapply(variance[["10X"]],
                                                       function(x) mean(x[intersect(names(x),
                                                                                    core.complex.list[[c]])])),
                                                value.name = "Var"),
                                           Method = "10X"),
                          cbind.data.frame(melt(lapply(variance[["Smart-seq2"]],
                                                       function(x) mean(x[intersect(names(x),
                                                                                    core.complex.list[[c]])])),
                                                value.name = "Var"),
                                           Method = "Smart-seq2"))
  plot.data <- rbind.data.frame(plot.data,
                                cbind.data.frame(merge(avg,var,"L1"),
                                                 "Complex" = c))
}

ggplot(plot.data) +
  geom_density(aes(x = Mean, color = Method.x)) +
  facet_wrap(~ L1)
ggplot(plot.data) +
  geom_density(aes(x = Var, color = Method.y)) +
  facet_wrap(~ L1)
# bc these are core complexes, avg expression is mid-high and var is low
# Smart-seq2 data looks better

pdf("Plots/0.5_Correlation/corecomplex_meanexpr_TS.pdf", width = 55, height = 25)
ggplot(plot.data) +
  geom_point(data = plot.data,
             aes(x = Complex, y = Mean, color = L1, shape = Method.x)) +
  scale_color_manual(values = cell_type_palette, name = "Cell Type") +
  xlab("Core complexes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 20), legend.position = "top")
dev.off()

# ------------------------------------------------------------------------------

# complex = "Histone H3.1 complex"
genes <- c("RPLP0","RPLP1")

plot.data <- data.frame()

for(file in files_SS2){
  
  cell_type <- gsub(".rds","",tail(strsplit(file,"_")[[1]],1))
  
  data <- readRDS(paste0("Outputs/3_TSDataPrep/smartseq2/Pseudobulk/Subsets/",
                         file))
  data <- log2(data+1)
  
  plot.data <- rbind.data.frame(plot.data, data.frame(t(data[genes,]),
                                                      "CellType" = cell_type))
}

ggplot(plot.data) +
  geom_point(aes_string(x = genes[1], y = genes[2],
                        color = "CellType"))



pvals <- c()

for(tissue in unique(plot.data$L1)){
  
  message(tissue)
  
  # Real p-value distribution
  tmp.foreground.data <- subset(foreground.data, (L1 == tissue) &
                                  (Var2 == "pval"))
  observed <- tmp.foreground.data$value
  
  # Permuted p-value distribution
  tmp.plot.data <- subset(plot.data, L1 == tissue)
  permuted <- tmp.plot.data$pval
  permuted <- sample(permuted, length(observed))
  
  # print(ks.test(x = observed, y = punif, alternative = "greater"))
  pvals <- c(pvals, ks.test(x = observed, y = permuted,
                            alternative = "greater")$p.value)
}

names(pvals) <- unique(plot.data$L1)
pvals <- p.adjust(pvals)
pvals <- sort(pvals)

# KS H0: two vectors come from the same distribution.
library(readr)

coreComplexes <- read_delim("Data/CORUM/coreComplexes.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
nrow(coreComplexes) # 3512

# Subset to human
coreComplexes <- subset(coreComplexes, Organism == "Human")
nrow(coreComplexes) # 2417

# Extract gene names per complex
core.complex.list <- c()
for(i in 1:nrow(coreComplexes)){
  core.complex.list[[coreComplexes$ComplexName[i]]] <- strsplit(coreComplexes$`subunits(Gene name)`[i],
                                                                split = ";")[[1]]
}
core.complex.list <- core.complex.list[lapply(core.complex.list, length) > 5] # 346


# Gene families
genefamilies <- read_delim("Data/GeneFamily/genefamilies.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)
genefamilies <- subset(genefamilies, !is.na(`Gene group name`))
genefamilies <- subset(genefamilies, !is.na(`Approved symbol`))
families <- unique(as.character(genefamilies$`Gene group name`)) # 3619
gene.family.list <- sapply(families,
                           function(fam) as.character(genefamilies$`Approved symbol`[genefamilies$`Gene group name` == fam]))
gene.family.list <- gene.family.list[unlist(lapply(gene.family.list, length)) > 5] # 612


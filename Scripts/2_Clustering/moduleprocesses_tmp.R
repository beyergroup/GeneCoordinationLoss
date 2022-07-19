# Presence of genes belonging to metabolic GO terms in biggest module

args <- commandArgs(trailingOnly = TRUE)
NET = args[1]
WDIR = args[2]
PATTERN = args[3]

.libPaths("Resources/Rlibs/R-4.0.3/")
library(biomaRt)
library(pheatmap)
source("Scripts/functions.R")


# load module calling
membership <- ReadRDS(paste0("Outputs/",WDIR,"/",NET,
    "_Adjacency_undirected_weightssum_rownormalized_membership.rds"))
modules <- names(table(membership)[table(membership) >= 10])

# focus on largest module
module = names(sort(table(membership), decreasing = T))[1]

# select GO terms of interest
GO_IDs <- c("cell differentiation" = "GO:0030154", # different linneages involved
            "cytoplasmic translation" = "GO:0002181", # components
            "signaling" = "GO:0023052", # signal transduction, signalling regulators, receptors
            "transcription, DNA-templated" = "GO:0006351", # extended transcription machinery and regulators
            
            "autophagy" = "GO:0006914", # members and regulators
            "carbohydrate metabolic process" = "GO:0005975", # pathway members
            "generation of precursor metabolites and energy" = "GO:0006091", # pathway members
            "ribosome biogenesis" = "GO:0042254", # includes ribosomal constituents and regulators
            
            "lysosome" = "GO:0005764", # lysosome components
            "mitochondrion" = "GO:0005739", # mito components
            "ribosome" = "GO:0005840", # ribosomal components
            
            "immune system process" = "GO:0002376", # e.g. immune cell activation
            "muscle system process" = "GO:0003012", # e.g. muscle contraction
            "nervous system process" = " GO:0050877", # e.g. synaptic transmission, members of GABA receptors
            "renal system process" = "GO:0003014", # e.g. filtration
            "reproductive process" = "GO:0022414", # e.g. meiosis
            
            # "cytoskeleton organization" = "GO:0007010", # assembly, disassembly, rearrangement of components
            # "extracellular matrix organization" = "GO:0030198", # assembly, disassembly, rearrangement of components
            # "lysosome organization" = "GO:0007040", # assembly, disassembly, rearrangement of components
            # "mitochondrion organization" = "GO:0007005", # biosynthesis, re-distribution, replication, regulation of mito
            # "peroxisome organization" = "GO:0007031", # assembly, disassembly, rearrangement of components
            NULL)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c("hgnc_symbol","go_id"),
               filters = "go", values = GO_IDs, mart = ensembl)
genes <- subset(genes, go_id %in% GO_IDs)
genes$GOTerm <- names(GO_IDs)[match(genes$go_id,GO_IDs)]

# count occurences across modules
genes <- subset(genes, hgnc_symbol %in% intersect(names(membership),genes$hgnc_symbol))
GO_IDs <- GO_IDs[GO_IDs %in% genes$go_id]

counts <- t(sapply(modules, function(m){
  members <- names(which(membership == m))
  sapply(names(GO_IDs), function(term) length(intersect(members,
                            subset(genes, GOTerm == term)$hgnc_symbol)))
}))
counts <- rbind("All" = sapply(names(GO_IDs),
                               function(term) length(intersect(names(membership),
                                                               subset(genes, GOTerm == term)$hgnc_symbol))),
                counts)

pheatmap(head(counts,21), scale = "none", legend = F,
         color = colorRampPalette(c("white","red"))(100),
         cluster_rows = F, cluster_cols = F,
         cellheight = 15, cellwidth = 20,
         display_numbers = head(counts,21),
         gaps_row = 1, gaps_col = c(8,11),
         filename = paste0("Plots/",WDIR,"/",NET,
              "_Adjacency_undirected_weightssum_rownormalized_GOslim.pdf"))

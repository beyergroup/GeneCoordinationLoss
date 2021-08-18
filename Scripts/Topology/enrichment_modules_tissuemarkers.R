# Enrichment in tissue-specific genes among module members

source("Scripts/functions.R")
library(limma, lib.loc = "Resources/Rlibs/R-4.0.3/")

NET = "stabsel_pcclasso"
METHOD = "greedy"
MODULE_FILE = paste0("Outputs/Human_Network/",NET,
                     "/Topology/Modules/membership_absweights_",
                     METHOD,"_iterreclustering.rds")


# read in network modules
modules <- readRDS(MODULE_FILE)

# read in tissue DE results
tissue_DE <- sapply(list.files("Outputs/GTEx/TissueDE", pattern = "sampled_tissueDE.rds", full.names = T),
                    readRDS, simplify = F)
names(tissue_DE) <- sapply(names(tissue_DE), function(x) strsplit(tail(strsplit(x, "/")[[1]],1), "_")[[1]][1])
tissue_DE <- lapply(tissue_DE, function(f) topTable(f, coef = 2, number = nrow(f$coefficients)))
tissue_DE <- lapply(tissue_DE, function(df) DataENSGToSymbol(df, remove_dup = T))

enrichments <- list()

# test for enrichment one module at the time
for(tissue in names(tissue_DE)){
  
  # GTEx background
  background <- rownames(tissue_DE[[tissue]])
  # tissue markers
  tissue_markers <- rownames(subset(tissue_DE[[tissue]], (logFC > 1) & (adj.P.Val < 0.05)))
  
  # net background
  background <- intersect(background, names(modules))
  tissue_markers <- intersect(tissue_markers, names(modules)) # tissue markers
  background <- setdiff(background, tissue_markers) # not tissue markers
  
  enrichments[[tissue]] <- data.frame()
  
  # iterate over net modules
  for(module in unique(modules)){
    module_genes <- names(which(modules == module))
    module_genes <- intersect(module_genes, union(tissue_markers,background))
    PP <- length(intersect(module_genes, tissue_markers)) # tissue markers AND module members
    NP <- length(intersect(module_genes, background)) # not tissue markers, but module members
    PN <- length(setdiff(tissue_markers, module_genes)) # tissue markers, but not module members
    NN <- length(setdiff(background,
                         union(tissue_markers,module_genes))) # neither tissue markers, nor module members
    
    if(length(tissue_markers) != (PP+PN))
      stop("Error in contingency table")
    
    if(length(module_genes) != (PP+NP))
      stop("Error in contingency table")
    
    if(length(intersect(background,names(modules))) != (NP+NN))
      stop("Error in contingency table")
    
    # create contingency table
    contingency <- matrix(c(PP,NP,PN,NN), nrow = 2,
                          dimnames = list("Tissue" = c("TissueMarker","NotTissueMarker"),
                                          "Module" = c("InsideModule","OutsideModule")))
    
    # test for independence
    Xsq <- chisq.test(contingency)
    # save p.value
    enrichments[[tissue]] <- rbind.data.frame(enrichments[[tissue]],
                                              data.frame("ModuleID" = module,
                                                         "pval" = Xsq$p.value))
  }
  
  # apply multiple test correction to p-values (#tests = #modules)
  enrichments[[tissue]]$padjust <- p.adjust(enrichments[[tissue]]$pval,
                                            method = "fdr")
}

saveRDS(enrichments, paste0("Outputs/Human_Network/",NET,
  "/Topology/Modules/TissueDE/enrichment_tissuemarkers_",METHOD,".rds"))

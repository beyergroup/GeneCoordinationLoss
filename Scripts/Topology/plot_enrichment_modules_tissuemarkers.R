# Explore further tissue-specific network modules

source("Scripts/functions.R")
library(reshape2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggplot2, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(igraph, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(ggpubr, lib.loc = "Resources/Rlibs/R-4.0.3/")
library(gridGraphics, lib.loc = "Resources/Rlibs/R-4.0.3/")

NET = "stabsel"
METHOD = "greedy"
MODULE_COL = "red"

enrichments <- readRDS(paste0("Outputs/Human_Network/",NET,
  "/Topology/Modules/TissueDE/enrichment_tissuemarkers_",METHOD,".rds"))

enrichments <- lapply(enrichments, function(df) subset(df, padjust < 0.05))
lapply(enrichments,dim) # less than 100
enrichments <- lapply(enrichments, function(df) df[order(df$padjust),])

# to what extent are these modules shared between tissues? Hopefully not much.
sort(table(unlist(lapply(enrichments, function(df) df$ModuleID))), decreasing = T)
# some are shared between up to 5 tissues - look at the enriched tissues + cross tissue
e <- melt(lapply(enrichments, function(df) df$ModuleID))

# read in membership info
membership <- readRDS(paste0("Outputs/Human_Network/",NET,
  "/Topology/Modules/membership_absweights_",METHOD,"_iterreclustering.rds"))

# read in GTEx expression levels
gtex_files <- list.files("GTEx_Networks/Tissue_Networks/Outputs", pattern = "sampled_data.rds",
                         full.names = T)
gtex_expr <- sapply(gtex_files, readRDS, simplify = F)
names(gtex_expr) <- sapply(names(gtex_expr),
  function(x) strsplit(tail(strsplit(x, "/")[[1]],1),"_")[[1]][1])
gtex_expr <- lapply(gtex_expr, DataENSGToSymbol)

# read in GO analysis
load(paste0("Outputs/Human_Network/",NET,
            "/Topology/Modules/GO_allcommunities_iterreclustering_",METHOD,".RData"))

# read in net itself
net_graph <- readRDS(paste0("Outputs/Human_Network/",NET,"/network_igraph.rds"))


for(tissue in names(enrichments)){
  
  modules <- head(enrichments[[tissue]],3)$ModuleID
  
  for(module in modules){
    
    plots <- list()
    
    genes <- names(membership)[which(membership == module)]
    tissues <- subset(e, value == module)$L1
    
    # Expression of module members across tissues
    plot.data <- melt(lapply(gtex_expr,
      function(m) rowMeans(m[intersect(genes,rownames(m)),])))
    plot.data$L1 <- factor(as.character(plot.data$L1),
      levels = c(unique(plot.data$L1)[-4],"CrossTissue"))
    plot.data$Enriched <- plot.data$L1 %in% tissues
    plots[[1]] <- ggplot(plot.data) +
      geom_boxplot(aes(x = L1, y = value), fill = "transparent") +
      geom_jitter(aes(x = L1, y = value, color = Enriched)) +
      scale_color_manual(values = c("TRUE" = MODULE_COL, "FALSE" = "black")) +
      ylab("Corrected expression") + guides(color = F) + xlab("") +
      ggtitle("Module expression across tissues") +
      theme(text = element_text(size = 20),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
    
    # GO enrichment
    plots[[2]] <- PlotGOEnrich(GOmf.list[[module]], col = MODULE_COL,
                               title = "GO molecular functions")
    
    pdf(paste0("Plots/Human_Network/",NET,"/Topology/Modules/",METHOD,"_",module,".pdf"),
        height = 9, width = 16)
    print(ggarrange(plotlist = plots, nrow = 1, widths = c(3,5), align = "h"))
    dev.off()
    
    # network module itself
    current_graph <- induced_subgraph(net_graph,
      vids = which(V(net_graph)$name %in% genes))
    E(current_graph)$width <- abs(E(current_graph)$weight)*2
    V(current_graph)$color <- MODULE_COL
    pdf(paste0("Plots/Human_Network/",NET,"/Topology/Modules/",METHOD,"_",module,"_netviz.pdf"))
    PlotNet(current_graph, title = paste("Module",module,METHOD), layout = "lgl")
    dev.off()
  }
}

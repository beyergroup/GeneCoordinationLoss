# Topology parameters

# Thresholds ------------------------------------------------------------------


# Color schemes ---------------------------------------------------------------

COLOR_SCHEME <- c("All genes" = "grey",
                  "Hub genes" = "darkviolet",
                  "Bottlenecks" = "darkmagenta",
                  "TFs" = "mediumvioletred",
                  "Poorly predicted" = "tomato2")

# File sizes ------------------------------------------------------------------

DIMS <- matrix(nrow = 2, ncol = 2, dimnames = list(c("All","Poorly Predicted"),
                                                   c("Height","Width")))
DIMS[1,] <- c(10,12)
DIMS[2,] <- c(12,12)


GO_PDF_DIMS <- expand.grid("Net" = c("stabsel","stabsel_filtered","stabsel_filtered_largestCC","stabsel_pcclasso",
                                     "stabsel_pcclasso_filtered","stabsel_pcclasso_filtered_largestCC"),
                           "Group" = c("Hub genes","Bottlenecks","Poorly predicted"))
GO_PDF_DIMS <- cbind(GO_PDF_DIMS,
                     "Height" = c(8,8,8,12,12,8,10,10,12,7,7,7,7,7,7,7,7,7),
                     "Width" = c(13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13))

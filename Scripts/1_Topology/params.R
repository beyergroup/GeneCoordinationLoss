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


GO_PDF_DIMS <- expand.grid("Net" = c("stabsel_filtered_trans_largestCC","stabsel_filtered_trans_largestCC_randomized",
                                     "stabsel_pcclasso_filtered_trans_largestCC", "stabsel_pcclasso_filtered_trans_largestCC_randomized"),
                           "Group" = c("Hub genes","Bottlenecks","Poorly predicted"))
GO_PDF_DIMS <- cbind(GO_PDF_DIMS,
                     "Height" = c(11,10,8,8,
                                  9,12,9,9,
                                  8,8,8,8),
                     "Width" = c(12,13,13,13,
                                 11,12,13,13,
                                 13,13,13,13))

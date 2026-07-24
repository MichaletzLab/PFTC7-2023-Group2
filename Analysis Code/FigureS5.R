# Figure S5: temperature sensitivities and reference values for E, gsw, and WUE
# vs PC2. Analogue of Fig. 5, generated from the same code but using PC2 rather 
# than PC1. All modelling and plotting code is in Figure5.R.
#
# Requires: DataPrep.R and PCA.R.

PC_AXIS <- "PC2"
FIG_TAG <- "FigureS5"

source("Analysis Code/Figure5.R")

message("\nPC_AXIS is now '", PC_AXIS, "'. Restart R, or set PC_AXIS <- \"PC1\" ",
        "and FIG_TAG <- \"Figure5\", before running Figure5.R in this session.")
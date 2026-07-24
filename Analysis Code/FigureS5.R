# Figure S5: temperature sensitivities and reference values for E, gsw, and WUE
# vs PC2. Analogue of Fig. 5, generated from the same code but using PC2 rather 
# than PC1. All modelling and plotting code is in Figure5.R.
#
# Requires: DataPrep.R and PCA.R.

PC_AXIS <- "PC2"

source("Analysis Code/Figure5.R")

message("\nPC_AXIS is now '", PC_AXIS, "'. Restart R, or set PC_AXIS <- \"PC1\", ",
        "before running the main-text figures in this session.")
# Figure S4: Sharpe-Schoolfield parameters for A and iWUE vs PC2.
# Analogue of Fig. 4, generated from same code but using PC2 rather than PC1.
# All modelling and plotting code is in Figure4.R.
#
# Requires: DataPrep.R, PCA.R, and outputs/figure4_ss_parameters.csv.

PC_AXIS <- "PC2"

source("Analysis Code/Figure4.R")

message("\nPC_AXIS is now '", PC_AXIS, "'. Restart R, or set PC_AXIS <- \"PC1\", ",
        "before running the main-text figures in this session.")
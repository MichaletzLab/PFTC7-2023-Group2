# ---- Fit GAMs with just Tleaf + PC1 + Species ----
gam_mod_A_PC1   <- gam(A    ~ s(Tleaf, k=3) + s(PC1, k=3) + ti(Tleaf, PC1, k=3) + Species,
                       data = raw.env.data_pca, method = "REML")
gam_mod_E_PC1   <- gam(E    ~ s(Tleaf, k=3) + s(PC1, k=3) + ti(Tleaf, PC1, k=3) + Species,
                       data = raw.env.data_pca, method = "REML")
gam_mod_gsw_PC1 <- gam(gsw  ~ s(Tleaf, k=3) + s(PC1, k=3) + ti(Tleaf, PC1, k=3) + Species,
                       data = raw.env.data_pca, method = "REML")
gam_mod_iWUE_PC1<- gam(iWUE ~ s(Tleaf, k=3) + s(PC1, k=3) + ti(Tleaf, PC1, k=3) + Species,
                       data = raw.env.data_pca, method = "REML")
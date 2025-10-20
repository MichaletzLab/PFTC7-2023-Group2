# ---- Fit GAMs ----
gam_mod_A_PC1   <- gam(A    ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5) 
                       + Species + s(curveID, bs="re"),                       
                       data = raw.env.data_pca, method = "REML")
#gam.check(gam_mod_A_PC1)

gam_mod_E_PC1   <- gam(E    ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5) + 
                         Species+ s(curveID, bs="re"),                       
                       data = raw.env.data_pca, method = "REML")
gam.check(gam_mod_E_PC1)

gam_mod_gsw_PC1 <- gam(gsw  ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5) + 
                         Species + s(curveID, bs="re"),                       
                       data = raw.env.data_pca, method = "REML")
#gam.check(gam_mod_gsw_PC1)

gam_mod_iWUE_PC1<- gam(iWUE ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5) + 
                         Species+ s(curveID, bs="re"),                       
                       data = raw.env.data_pca, method = "REML")
#gam.check(gam_mod_iWUE_PC1)

gam_mod_eWUE_PC1<- gam(eWUE ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5) + 
                         Species+ s(curveID, bs="re"),                       
                       data = raw.env.data_pca, method = "REML")



summary(gam_mod_A_PC1)
summary(gam_mod_E_PC1)
summary(gam_mod_gsw_PC1)
summary(gam_mod_iWUE_PC1)
summary(gam_mod_eWUE_PC1)



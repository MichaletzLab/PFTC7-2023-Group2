# libraries
library(mgcv)

# For Fig 3 GAMs, PC1 is a site-level predictor, it takes one value per site
# (9 total), so curves and points within a site aren't independent. Without a 
# site-level term, the tests on s(PC1) and ti(Tleaf, PC1) treat within-site 
# observatins as independent (pseudoreplication). s(Tleaf, Site, bs="fs")
# gives each site its own temperature response, taking care of site-level 
# variation so s(PC1) and ti(Tleaf, PC1) are tested against real between-site
# variation. This is similar to how site random effects added in Fig. 4 and 5
# mixed models.

# Site corresponds to unique elevations, which is grouping that PC1 varies at
# (for the factor smooth)
raw.env.data_pca$Site <- factor(raw.env.data_pca$Elevation)

gam_mod_A_PC1 <- gam(A ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5)
                     + Species + s(curveID, bs="re") + s(Tleaf, Site, bs="fs", k=5),
                     data = raw.env.data_pca, method = "REML")

gam_mod_E_PC1 <- gam(E ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5)
                     + Species + s(curveID, bs="re") + s(Tleaf, Site, bs="fs", k=5),
                     data = raw.env.data_pca, method = "REML")

gam_mod_gsw_PC1 <- gam(gsw ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5)
                       + Species + s(curveID, bs="re") + s(Tleaf, Site, bs="fs", k=5),
                       data = raw.env.data_pca, method = "REML")

gam_mod_iWUE_PC1 <- gam(iWUE ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5)
                        + Species + s(curveID, bs="re") + s(Tleaf, Site, bs="fs", k=5),
                        data = raw.env.data_pca, method = "REML")

gam_mod_WUE_PC1 <- gam(WUE ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5)
                       + Species + s(curveID, bs="re") + s(Tleaf, Site, bs="fs", k=5),
                       data = raw.env.data_pca, method = "REML")

# ---- Fit GAMs ----
gam_mod_A_PC1   <- gam(A    ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5) 
                       + Species + s(curveID, bs="re"),                       
                       data = raw.env.data_pca, method = "REML")
summary(gam_mod_A_PC1)
#gam.check(gam_mod_A_PC1)

gam_mod_E_PC1   <- gam(E    ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5) + 
                         Species+ s(curveID, bs="re"),                       
                       data = raw.env.data_pca, method = "REML")
#gam.check(gam_mod_E_PC1)

gam_mod_gsw_PC1 <- gam(gsw  ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5) + 
                         Species + s(curveID, bs="re"),                       
                       data = raw.env.data_pca, method = "REML")
#gam.check(gam_mod_gsw_PC1)

gam_mod_iWUE_PC1<- gam(iWUE ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5) + 
                         Species+ s(curveID, bs="re"),                       
                       data = raw.env.data_pca, method = "REML")
#gam.check(gam_mod_iWUE_PC1)

gam_mod_WUE_PC1<- gam(WUE ~ s(Tleaf, k=5) + s(PC1, k=5) + ti(Tleaf, PC1, k=5) + 
                         Species+ s(curveID, bs="re"),                       
                       data = raw.env.data_pca, method = "REML")
#gam.check(gam_mod_WUE_PC1)

summary(gam_mod_A_PC1)
summary(gam_mod_E_PC1)
summary(gam_mod_gsw_PC1)
summary(gam_mod_iWUE_PC1)
summary(gam_mod_WUE_PC1)


#For supplemental
gam_mod_full_environ_A <- gam(A ~ s(Tleaf, k=4) + 
                               s(mean_moist_pct, k=4) + 
                               s(mean_T1, k=4) + 
                               s(mean_T2, k=4) + 
                               s(mean_air, k=4) + 
                               s(vegetation_height, k=4) + 
                        Species+ s(curveID, bs="re"),                       
                      data = raw.env.data_pca, method = "REML")
summary(gam_mod_full_environ_A)
gam_mod_full_environ_E <- gam(E ~ s(Tleaf, k=4) + 
                               s(mean_moist_pct, k=4) + 
                               s(mean_T1, k=4) + 
                               s(mean_T2, k=4) + 
                               s(mean_air, k=4) + 
                               s(vegetation_height, k=4) + 
                               Species+ s(curveID, bs="re"),                       
                             data = raw.env.data_pca, method = "REML")
gam_mod_full_environ_gsw <- gam(gsw ~ s(Tleaf, k=4) + 
                               s(mean_moist_pct, k=4) +
                               s(mean_T1, k=4) +
                               s(mean_T2, k=4) + 
                               s(mean_air, k=4) + 
                               s(vegetation_height, k=4) + 
                               Species+ s(curveID, bs="re"),                       
                             data = raw.env.data_pca, method = "REML")
gam_mod_full_environ_iWUE <- gam(iWUE ~ s(Tleaf, k=4) + 
                                 s(mean_moist_pct, k=4) +
                                 s(mean_T1, k=4) + 
                                 s(mean_T2, k=4) + 
                                 s(mean_air, k=4) + 
                                 s(vegetation_height, k=4) + 
                                 Species+ s(curveID, bs="re"),                       
                               data = raw.env.data_pca, method = "REML")
gam_mod_full_environ_WUE <- gam(WUE ~ s(Tleaf, k=4) + 
                                   s(mean_moist_pct, k=4) + 
                                   s(mean_T1, k=4) + 
                                   s(mean_T2, k=4) + 
                                   s(mean_air, k=4) +  
                                   s(vegetation_height, k=4) + 
                                   Species+ s(curveID, bs="re"),                       
                                 data = raw.env.data_pca, method = "REML")
gam_mod_full_PC_A <- gam(A ~ s(Tleaf, k=4) + 
                                   s(PC1, k=4) +
                                   s(PC2, k=4) +
                                   s(PC3, k=4) + 
                                   s(PC4, k=4) + 
                                   s(PC5, k=4) +  
                                   Species+ s(curveID, bs="re"),                       
                                 data = raw.env.data_pca, method = "REML")
gam_mod_full_PC_E <- gam(E ~ s(Tleaf, k=4) + 
                           s(PC1, k=4) + 
                           s(PC2, k=4) + 
                           s(PC3, k=4) + 
                           s(PC4, k=4) + 
                           s(PC5, k=4) + 
                           Species+ s(curveID, bs="re"),                       
                         data = raw.env.data_pca, method = "REML")
gam_mod_full_PC_gsw <- gam(gsw ~ s(Tleaf, k=4) + 
                           s(PC1, k=4) + 
                           s(PC2, k=4) + 
                           s(PC3, k=4) +
                           s(PC4, k=4) + 
                           s(PC5, k=4) +  
                           Species+ s(curveID, bs="re"),                       
                         data = raw.env.data_pca, method = "REML")
gam_mod_full_PC_iWUE <- gam(iWUE ~ s(Tleaf, k=4) + 
                           s(PC1, k=4) + 
                           s(PC2, k=4) +
                           s(PC3, k=4) +
                           s(PC4, k=4) + 
                           s(PC5, k=4) +  
                           Species+ s(curveID, bs="re"),                       
                         data = raw.env.data_pca, method = "REML")
gam_mod_full_PC_WUE <- gam(WUE ~ s(Tleaf, k=4) + 
                              s(PC1, k=4) + 
                              s(PC2, k=4) +
                              s(PC3, k=4) +  
                              s(PC4, k=4) + 
                              s(PC5, k=4) + 
                              Species+ s(curveID, bs="re"),                       
                            data = raw.env.data_pca, method = "REML")




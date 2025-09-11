raw.env.data_pca$Country <- factor(raw.env.data_pca$Country)
library(mgcv)
raw.env.data_pca$Country <- factor(raw.env.data_pca$Country)

# All PCs + Elevation
gam_mod_all_Elev <- gam(A ~ s(Tleaf, k=3) +
                          s(Elevation, k=3) +
                          s(PC1, k=3) + s(PC2, k=3) + s(PC3, k=3) +
                          ti(Tleaf, Elevation, k=3) +
                          ti(Tleaf, PC1, k=3) + ti(Tleaf, PC2, k=3) +
                          ti(Tleaf, PC3, k=3),
                        data = raw.env.data_pca,
                        method="REML")

# ️PC1 + Elevation
gam_mod_PC1_Elev <- gam(A ~ s(Tleaf, k=3) +
                          s(Elevation, k=3) +
                          s(PC1, k=3) +
                          ti(Tleaf, Elevation, k=3) +
                          ti(Tleaf, PC1, k=3),
                        data = raw.env.data_pca,
                        method="REML")

# All PCs + Country
gam_mod_all_Country <- gam(A ~ s(Tleaf, k=3) +
                             Country +
                             s(PC1, k=3) + s(PC2, k=3) + s(PC3, k=3) +
                             ti(Tleaf, PC1, k=3) + ti(Tleaf, PC2, k=3) +
                             ti(Tleaf, PC3, k=3),
                           data = raw.env.data_pca,
                           method="REML")

# PC1 + Country
gam_mod_PC1_Country <- gam(A ~ s(Tleaf, k=3) +
                             Country +
                             s(PC1, k=3) +
                             ti(Tleaf, PC1, k=3),
                           data = raw.env.data_pca,
                           method="REML")

# Summaries
summary(gam_mod_all_Elev)
summary(gam_mod_PC1_Elev)
summary(gam_mod_all_Country)
summary(gam_mod_PC1_Country)


library(ggplot2)

plot_gam_response <- function(model, response="A", include="Elevation") {
  # create newdata for predictions
  newdat <- with(raw.env.data_pca,
                 data.frame(
                   Tleaf = seq(min(Tleaf, na.rm=TRUE),
                               max(Tleaf, na.rm=TRUE),
                               length.out=200),
                   PC1 = mean(PC1, na.rm=TRUE),
                   PC2 = mean(PC2, na.rm=TRUE),
                   PC3 = mean(PC3, na.rm=TRUE),
                   PC4 = mean(PC4, na.rm=TRUE),
                   Elevation = median(Elevation, na.rm=TRUE),
                   Country = factor("Norway", levels=levels(raw.env.data_pca$Country))
                 ))
  
  # predictions
  p <- predict(model, newdata=newdat, se.fit=TRUE)
  newdat$fit <- p$fit
  newdat$se <- p$se.fit
  
  ggplot(newdat, aes(x=Tleaf, y=fit)) +
    geom_point(data=raw.env.data_pca,
               aes(x=Tleaf, y=.data[[response]]),
               inherit.aes=FALSE, alpha=0.01, size=0.5) +
    geom_line(color="darkgreen", linewidth=1) +
    geom_ribbon(aes(ymin=fit-2*se, ymax=fit+2*se),
                fill="lightgreen", alpha=0.3) +
    theme_classic() +
    labs(x="Tleaf", y=response,
         title=paste("GAM fit for", response, "(", include, ")"))
}
p_all_Elev     <- plot_gam_response(gam_mod_all_Elev, "A", include="Elevation")
p_PC1_Elev     <- plot_gam_response(gam_mod_PC1_Elev, "A", include="Elevation")
p_all_Country  <- plot_gam_response(gam_mod_all_Country, "A", include="Country")
p_PC1_Country  <- plot_gam_response(gam_mod_PC1_Country, "A", include="Country")

p_all_Elev | p_PC1_Elev | p_all_Country | p_PC1_Country


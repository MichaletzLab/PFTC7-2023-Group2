# Check correlation between Elevation and T2 (note using T2 becuase of PCA results)----
cor.test(raw.env.data$Elevation, raw.env.data$mean_T2, method = "pearson") #This is not very correlated! If you replace this with the summary of elevation that just has the 9 points and 9 T2 values, there is no significant correlation.
cor.test(raw.env.data$Elevation, raw.env.data$mean_moist_pct, method = "pearson")
# Full A GAM -------------------------------------------------------------------
gam_mod.A_species <- gam(A ~ s(Tleaf, k = 3) +
                           s(Elevation, k = 3) +
                           s(mean_T2, k = 3) +
                           s(mean_moist_pct, k = 3) +
                           Species +
                           ti(Tleaf, Elevation, k = 3) +
                           ti(Tleaf, mean_T2, k = 3) +
                           ti(Tleaf, mean_moist_pct, k = 3) +
                           s(Tleaf, by = Species, k = 3),
                         data = raw.env.data,
                         method = "REML"
)
summary(gam_mod.A_species)
# Plot A GAM --------------------------------------------------------------------
gam_mod.A <- gam(A ~ s(Tleaf, k = 3) +
                   s(Elevation, k = 3) +
                   s(mean_T2, k = 3) +
                   s(mean_moist_pct, k = 3) +
                   ti(Tleaf, Elevation, k = 3) +
                   ti(Tleaf, mean_T2, k = 3) +
                   ti(Tleaf, mean_moist_pct, k = 3),
                 data = raw.env.data,
                 method = "REML"
)
summary(gam_mod.A)

make_preds <- function(model, data, driver, response) {
  qs <- quantile(data[[driver]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
  names(qs) <- c("Low", "Medium", "High")
  
  out <- lapply(names(qs), function(level) {
    newdat <- expand.grid(
      Tleaf = seq(min(data$Tleaf, na.rm=TRUE),
                  max(data$Tleaf, na.rm=TRUE),
                  length.out=200),
      Elevation = median(data$Elevation, na.rm=TRUE),
      mean_T2 = median(data$mean_T2, na.rm=TRUE),
      mean_moist_pct = median(data$mean_moist_pct, na.rm=TRUE)
    )
    newdat[[driver]] <- qs[level]
    
    p <- predict(model, newdata=newdat, se.fit=TRUE)
    newdat$fit <- p$fit
    newdat$se <- p$se.fit
    newdat$driver <- driver
    newdat$level <- level
    newdat$response <- response
    newdat
  })
  
  bind_rows(out)
}
preds_elev  <- make_preds(gam_mod.A, raw.env.data, "Elevation", "A")
preds_T2    <- make_preds(gam_mod.A, raw.env.data, "mean_T2", "A")
preds_moist <- make_preds(gam_mod.A, raw.env.data, "mean_moist_pct", "A")

plot_driver <- function(preds, driver_label, data, driver) {
  
  # grab the response name from preds (all rows have the same value)
  response <- unique(preds$response)
  
  # compute cutpoints (same as in make_preds)
  qs <- quantile(data[[driver]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
  names(qs) <- c("Low", "Medium", "High")
  
  # assign each row in raw data to Low/Medium/High
  data$level <- cut(data[[driver]],
                    breaks = c(-Inf, qs["Low"], qs["Medium"], qs["High"], Inf),
                    labels = c("Low", "Medium", "High", NA),
                    include.lowest = TRUE, right = TRUE)
  
  # make plot
  ggplot(preds, aes(x=Tleaf, y=fit, color=level)) +
    # faint raw points, colored by level
    geom_point(data = data, aes(x=Tleaf, y=.data[[response]], color=level),
               inherit.aes = FALSE, alpha = 0.005, size = 0.5) +
    geom_line(size=1) +
    geom_ribbon(aes(ymin=fit-2*se, ymax=fit+2*se, fill=level), 
                alpha=0.2, colour=NA) +
    scale_color_manual(values=c("Low"="lightgreen", "Medium"="green3", "High"="darkgreen"),
                       breaks=c("Low","Medium","High")) +
    scale_fill_manual(values=c("Low"="lightgreen", "Medium"="green3", "High"="darkgreen"),
                      breaks=c("Low","Medium","High")) +
    theme_classic() +
    labs(x="Tleaf", 
         y=response, 
         color=driver_label, 
         fill=driver_label, 
         title=NULL)
}



p1_A <- plot_driver(preds_elev, "Elevation", raw.env.data, "Elevation")
p2_A <- plot_driver(preds_T2, "Temperature", raw.env.data, "mean_T2")
p3_A <- plot_driver(preds_moist, "Moisture", raw.env.data, "mean_moist_pct")


(p1_A | p2_A | p3_A)


# Full E GAM -------------------------------------------------------------------
gam_mod.E_species <- gam(E ~ s(Tleaf, k = 3) +
                           s(Elevation, k = 3) +
                           s(mean_T2, k = 3) +
                           s(mean_moist_pct, k = 3) +
                           Species +
                           ti(Tleaf, Elevation, k = 3) +
                           ti(Tleaf, mean_T2, k = 3) +
                           ti(Tleaf, mean_moist_pct, k = 3) +
                           s(Tleaf, by = Species, k = 3),
                         data = raw.env.data,
                         method = "REML"
)
summary(gam_mod.E_species)

# Plot E GAM -------------------------------------------------------------------
gam_mod.E <- gam(E ~ s(Tleaf, k = 3) +
                   s(Elevation, k = 3) +
                   s(mean_T2, k = 3) +
                   s(mean_moist_pct, k = 3) +
                   ti(Tleaf, Elevation, k = 3) +
                   ti(Tleaf, mean_T2, k = 3) +
                   ti(Tleaf, mean_moist_pct, k = 3),
                 data = raw.env.data,
                 method = "REML"
)
summary(gam_mod.E)

# Predictions for E
preds_elev_E  <- make_preds(gam_mod.E, raw.env.data, "Elevation", "E")
preds_T2_E    <- make_preds(gam_mod.E, raw.env.data, "mean_T2", "E")
preds_moist_E <- make_preds(gam_mod.E, raw.env.data, "mean_moist_pct", "E")

# Plots for E

p1_E <- plot_driver(preds_elev_E, "Elevation", raw.env.data, "Elevation")
p2_E <- plot_driver(preds_T2_E, "Temperature", raw.env.data, "mean_T2")
p3_E <- plot_driver(preds_moist_E, "Moisture", raw.env.data, "mean_moist_pct")

(p1_E | p2_E | p3_E)



# Full gsw GAM -----------------------------------------------------------------
gam_mod.gsw_species <- gam(gsw ~ s(Tleaf, k = 3) +
                             s(Elevation, k = 3) +
                             s(mean_T2, k = 3) +
                             s(mean_moist_pct, k = 3) +
                             Species +
                             ti(Tleaf, Elevation, k = 3) +
                             ti(Tleaf, mean_T2, k = 3) +
                             ti(Tleaf, mean_moist_pct, k = 3) +
                             s(Tleaf, by = Species, k = 3),
                           data = raw.env.data,
                           method = "REML"
)
summary(gam_mod.gsw_species)

# Plot gsw GAM -----------------------------------------------------------------
gam_mod.gsw <- gam(gsw ~ s(Tleaf, k = 3) +
                     s(Elevation, k = 3) +
                     s(mean_T2, k = 3) +
                     s(mean_moist_pct, k = 3) +
                     ti(Tleaf, Elevation, k = 3) +
                     ti(Tleaf, mean_T2, k = 3) +
                     ti(Tleaf, mean_moist_pct, k = 3),
                   data = raw.env.data,
                   method = "REML"
)
summary(gam_mod.gsw)

# Predictions for gsw
preds_elev_gsw  <- make_preds(gam_mod.gsw, raw.env.data, "Elevation", "gsw")
preds_T2_gsw    <- make_preds(gam_mod.gsw, raw.env.data, "mean_T2", "gsw")
preds_moist_gsw <- make_preds(gam_mod.gsw, raw.env.data, "mean_moist_pct", "gsw")

# Plots for gsw
p1_gsw <- plot_driver(preds_elev_gsw, "Elevation", raw.env.data, "Elevation")
p2_gsw <- plot_driver(preds_T2_gsw, "Temperature", raw.env.data, "mean_T2")
p3_gsw <- plot_driver(preds_moist_gsw, "Moisture", raw.env.data, "mean_moist_pct")

(p1_gsw | p2_gsw | p3_gsw)

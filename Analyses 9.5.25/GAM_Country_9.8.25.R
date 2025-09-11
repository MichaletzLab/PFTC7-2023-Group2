#-------------------------------
# Full A GAM with Species
#-------------------------------
gam_mod.A_species_c <- gam(A ~ s(Tleaf, k = 3) +
                             Country +
                             s(mean_T2, k = 3) +
                             s(mean_moist_pct, k = 3) +
                             Species +
                             ti(Tleaf, mean_T2, k = 3) +
                             ti(Tleaf, mean_moist_pct, k = 3) +
                             s(Tleaf, by = Species, k = 3),
                           data = raw.env.data,
                           method = "REML")
summary(gam_mod.A_species_c) #Not sure I can have species and country in the same model

#-------------------------------
# Plot A GAM (no Species)
#-------------------------------
gam_mod.A_c <- gam(A ~ s(Tleaf, k = 3) +
                     Country +
                     s(mean_T2, k = 3) +
                     s(mean_moist_pct, k = 3) +
                     ti(Tleaf, mean_T2, k = 3) +
                     ti(Tleaf, mean_moist_pct, k = 3),
                   data = raw.env.data,
                   method = "REML")
summary(gam_mod.A_c)

#-------------------------------
# Prediction helper
#-------------------------------
raw.env.data$Country <- factor(raw.env.data$Country)

make_preds_c <- function(model, data, driver, response) {
  if (driver == "Country") {
    levels_driver <- levels(data$Country)
    out <- lapply(levels_driver, function(level) {
      newdat <- expand.grid(
        Tleaf = seq(min(data$Tleaf, na.rm=TRUE),
                    max(data$Tleaf, na.rm=TRUE),
                    length.out=200),
        Country = factor(level, levels = levels(data$Country)),
        mean_T2 = median(data$mean_T2, na.rm=TRUE),
        mean_moist_pct = median(data$mean_moist_pct, na.rm=TRUE)
      )
      p <- predict(model, newdata=newdat, se.fit=TRUE)
      newdat$fit <- p$fit
      newdat$se <- p$se.fit
      newdat$driver <- driver
      newdat$level <- level
      newdat$response <- response
      newdat
    })
  } else {
    qs <- quantile(data[[driver]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
    names(qs) <- c("Low", "Medium", "High")
    out <- lapply(names(qs), function(level) {
      newdat <- expand.grid(
        Tleaf = seq(min(data$Tleaf, na.rm=TRUE),
                    max(data$Tleaf, na.rm=TRUE),
                    length.out=200),
        Country = factor("Norway", levels = levels(data$Country)), # default
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
  }
  bind_rows(out)
}

#-------------------------------
# Plotting helper
#-------------------------------
library(rlang)

plot_driver_c <- function(preds, driver_label, data, driver) {
  response <- unique(preds$response)  # "A", "E", etc.
  
  if (driver == "Country") {
    ggplot(preds, aes(x=Tleaf, y=fit, color=level)) +
      geom_point(data=data, aes(x=Tleaf, y=!!sym(response), color=Country),
                 inherit.aes=FALSE, alpha=0.01, size=0.5) +
      geom_line(linewidth=1) +
      geom_ribbon(aes(ymin=fit-2*se, ymax=fit+2*se, fill=level),
                  alpha=0.2, colour=NA) +
      theme_classic() +
      labs(x="Tleaf", y=response, color=driver_label, fill=driver_label, title=NULL)
  } else {
    qs <- quantile(data[[driver]], probs=c(0.1,0.5,0.9), na.rm=TRUE)
    names(qs) <- c("Low","Medium","High")
    data$level <- cut(data[[driver]],
                      breaks=c(-Inf, qs["Low"], qs["Medium"], qs["High"], Inf),
                      labels=c("Low","Medium","High", NA),
                      include.lowest=TRUE, right=TRUE)
    
    ggplot(preds, aes(x=Tleaf, y=fit, color=level)) +
      geom_point(data=data, aes(x=Tleaf, y=!!sym(response), color=level),
                 inherit.aes=FALSE, alpha=0.01, size=0.5) +
      geom_line(linewidth=1) +
      geom_ribbon(aes(ymin=fit-2*se, ymax=fit+2*se, fill=level),
                  alpha=0.2, colour=NA) +
      scale_color_manual(values=c("Low"="lightgreen","Medium"="green3","High"="darkgreen"),
                         breaks=c("Low","Medium","High")) +
      scale_fill_manual(values=c("Low"="lightgreen","Medium"="green3","High"="darkgreen"),
                        breaks=c("Low","Medium","High")) +
      theme_classic() +
      labs(x="Tleaf", y=response, color=driver_label, fill=driver_label, title=NULL)
  }
}

#-------------------------------
# Predictions and plots for A
#-------------------------------
preds_country_A_c <- make_preds_c(gam_mod.A_c, raw.env.data, "Country", "A")
preds_T2_A_c      <- make_preds_c(gam_mod.A_c, raw.env.data, "mean_T2", "A")
preds_moist_A_c   <- make_preds_c(gam_mod.A_c, raw.env.data, "mean_moist_pct", "A")

p1_A_c <- plot_driver_c(preds_country_A_c, "Country", raw.env.data, "Country")
p2_A_c <- plot_driver_c(preds_T2_A_c, "Temperature", raw.env.data, "mean_T2")
p3_A_c <- plot_driver_c(preds_moist_A_c, "Moisture", raw.env.data, "mean_moist_pct")

(p1_A_c | p2_A_c | p3_A_c)

#-------------------------------
# Repeat for E
#-------------------------------
gam_mod.E_c <- gam(E ~ s(Tleaf, k = 3) +
                     Country +
                     s(mean_T2, k = 3) +
                     s(mean_moist_pct, k = 3) +
                     ti(Tleaf, mean_T2, k = 3) +
                     ti(Tleaf, mean_moist_pct, k = 3),
                   data = raw.env.data,
                   method = "REML")

preds_country_E_c <- make_preds_c(gam_mod.E_c, raw.env.data, "Country", "E")
preds_T2_E_c      <- make_preds_c(gam_mod.E_c, raw.env.data, "mean_T2", "E")
preds_moist_E_c   <- make_preds_c(gam_mod.E_c, raw.env.data, "mean_moist_pct", "E")

p1_E_c <- plot_driver_c(preds_country_E_c, "Country", raw.env.data, "Country")
p2_E_c <- plot_driver_c(preds_T2_E_c, "Temperature", raw.env.data, "mean_T2")
p3_E_c <- plot_driver_c(preds_moist_E_c, "Moisture", raw.env.data, "mean_moist_pct")

(p1_E_c | p2_E_c | p3_E_c)

#-------------------------------
# Repeat for gsw
#-------------------------------
gam_mod.gsw_c <- gam(gsw ~ s(Tleaf, k = 3) +
                       Country +
                       s(mean_T2, k = 3) +
                       s(mean_moist_pct, k = 3) +
                       ti(Tleaf, mean_T2, k = 3) +
                       ti(Tleaf, mean_moist_pct, k = 3),
                     data = raw.env.data,
                     method = "REML")

preds_country_gsw_c <- make_preds_c(gam_mod.gsw_c, raw.env.data, "Country", "gsw")
preds_T2_gsw_c      <- make_preds_c(gam_mod.gsw_c, raw.env.data, "mean_T2", "gsw")
preds_moist_gsw_c   <- make_preds_c(gam_mod.gsw_c, raw.env.data, "mean_moist_pct", "gsw")

p1_gsw_c <- plot_driver_c(preds_country_gsw_c, "Country", raw.env.data, "Country")
p2_gsw_c <- plot_driver_c(preds_T2_gsw_c, "Temperature", raw.env.data, "mean_T2")
p3_gsw_c <- plot_driver_c(preds_moist_gsw_c, "Moisture", raw.env.data, "mean_moist_pct")

(p1_gsw_c | p2_gsw_c | p3_gsw_c)

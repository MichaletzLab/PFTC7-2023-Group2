raw.env.data <- raw.env.data %>%
  mutate(iWUE = A/gsw)

# 2. Make categories
elev_q <- quantile(raw.env.data$Elevation, probs = c(0.33,0.66), na.rm=TRUE)
t2_q   <- quantile(raw.env.data$mean_T2, probs = c(0.33,0.66), na.rm=TRUE)
moist_q<- quantile(raw.env.data$mean_moist_pct, probs = c(0.33,0.66), na.rm=TRUE)

raw.env.data <- raw.env.data %>%
  mutate(
    Elev_cat = cut(Elevation, breaks = c(-Inf,elev_q,Inf), labels=c("Low","Medium","High")),
    T2_cat   = cut(mean_T2, breaks = c(-Inf,t2_q,Inf), labels=c("Low","Medium","High")),
    Moist_cat= cut(mean_moist_pct, breaks = c(-Inf,moist_q,Inf), labels=c("Low","Medium","High"))
  )
gam_iWUE <- gam(iWUE ~ s(Tleaf,k=3) +
                  s(Elevation,k=3) +
                  s(mean_T2,k=3) +
                  s(mean_moist_pct,k=3) +
                  ti(Tleaf,Elevation,k=3) +
                  ti(Tleaf,mean_T2,k=3) +
                  ti(Tleaf,mean_moist_pct,k=3),
                data = raw.env.data,
                method = "REML")
summary(gam_iWUE)

gam_iWUE_Spe <- gam(iWUE ~ s(Tleaf,k=3) +
                  s(Elevation,k=3) +
                  s(mean_T2,k=3) +
                  Species+
                  s(mean_moist_pct,k=3) +
                  ti(Tleaf,Elevation,k=3) +
                  ti(Tleaf,mean_T2,k=3) +
                  s(Tleaf, by = Species, k = 3)+
                  ti(Tleaf,mean_moist_pct,k=3),
                data = raw.env.data,
                method = "REML")
summary(gam_iWUE_Spe)

make_preds <- function(model, data, driver, response) {
  qs <- quantile(data[[driver]], probs = c(0.1,0.5,0.9), na.rm=TRUE)
  names(qs) <- c("Low","Medium","High")
  
  out <- lapply(names(qs), function(level){
    newdat <- expand.grid(
      Tleaf = seq(min(data$Tleaf,na.rm=TRUE),
                  max(data$Tleaf,na.rm=TRUE), length.out=200),
      Elevation = median(data$Elevation,na.rm=TRUE),
      mean_T2   = median(data$mean_T2,na.rm=TRUE),
      mean_moist_pct = median(data$mean_moist_pct,na.rm=TRUE)
    )
    newdat[[driver]] <- qs[level]
    p <- predict(model,newdata=newdat,se.fit=TRUE)
    newdat$fit <- p$fit; newdat$se <- p$se.fit
    newdat$driver <- driver; newdat$level <- level; newdat$response <- response
    newdat
  })
  bind_rows(out)
}

preds_elev_iWUE  <- make_preds(gam_iWUE, raw.env.data, "Elevation",      "iWUE")
preds_T2_iWUE    <- make_preds(gam_iWUE, raw.env.data, "mean_T2",        "iWUE")
preds_moist_iWUE <- make_preds(gam_iWUE, raw.env.data, "mean_moist_pct", "iWUE")

# Purple gradient colors
purple_cols <- c("Low"="#d3bdf0", "Medium"="#9b59b6", "High"="#5e3370")

# Plotting function
plot_driver <- function(preds, driver_label, data, driver) {
  response <- unique(preds$response)
  qs <- quantile(data[[driver]], probs = c(0.1,0.5,0.9), na.rm=TRUE)
  names(qs) <- c("Low","Medium","High")
  data$level <- cut(data[[driver]],
                    breaks=c(-Inf,qs["Low"],qs["Medium"],qs["High"],Inf),
                    labels=c("Low","Medium","High",NA),
                    include.lowest=TRUE,right=TRUE)
  
  ggplot(preds, aes(x=Tleaf, y=fit, color=level)) +
    geom_point(data=data, aes(x=Tleaf,y=.data[[response]],color=level),
               inherit.aes=FALSE, alpha=0.01, size=0.5) +
    geom_line(size=1) +
    geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=level),
                alpha=0.2,colour=NA) +
    scale_color_manual(values=purple_cols, breaks=c("Low","Medium","High")) +
    scale_fill_manual(values=purple_cols, breaks=c("Low","Medium","High")) +
    theme_classic() +
    coord_cartesian(ylim=c(0,250)) +
    labs(
      x="Tleaf [°C]",
      y="iWUE [µmol mol⁻¹]",
      color=driver_label,
      fill=driver_label
    )
}

# Three-panel figure
p1_iWUE <- plot_driver(preds_elev_iWUE,  "Elevation",   raw.env.data, "Elevation")
p2_iWUE <- plot_driver(preds_T2_iWUE,    "Temperature", raw.env.data, "mean_T2")
p3_iWUE <- plot_driver(preds_moist_iWUE, "Moisture",    raw.env.data, "mean_moist_pct")

ggarrange(p1_iWUE, p2_iWUE, p3_iWUE, nrow=1, labels=c("A","B","C"))

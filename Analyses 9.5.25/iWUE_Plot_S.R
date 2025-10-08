# ---- Fit GAM with Species as a factor ----
gam_iWUE_Spe <- gam(iWUE ~ 
                      s(Tleaf, k = 4) +
                      s(Elevation, k = 4) +
                      s(mean_T2, k = 4) +
                      s(mean_moist_pct, k = 4) +
                      ti(Tleaf, Elevation, k = 4) +
                      ti(Tleaf, mean_T2, k = 4) +
                      ti(Tleaf, mean_moist_pct, k = 4) +
                      Species +                          # parametric species effect
                      s(Tleaf, by = Species, k = 4),     # species-specific smooths
                    data = raw.env.data,
                    method = "REML"
)

# ---- Choose a reference species ----
target_species <- names(sort(table(raw.env.data$Species), decreasing = TRUE))[1]

# ---- Prediction function ----
make_preds_with_species <- function(model, data, driver, response, target_species) {
  qs <- quantile(data[[driver]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
  names(qs) <- c("Low", "Medium", "High")
  
  out <- lapply(names(qs), function(level) {
    newdat <- expand.grid(
      Tleaf = seq(min(data$Tleaf, na.rm = TRUE),
                  max(data$Tleaf, na.rm = TRUE), length.out = 200),
      Elevation = median(data$Elevation, na.rm = TRUE),
      mean_T2   = median(data$mean_T2, na.rm = TRUE),
      mean_moist_pct = median(data$mean_moist_pct, na.rm = TRUE),
      Species = factor(target_species, levels = levels(data$Species))
    )
    newdat[[driver]] <- qs[level]
    
    p <- predict(model, newdata = newdat, se.fit = TRUE)
    newdat$fit <- p$fit
    newdat$se <- p$se.fit
    newdat$driver <- driver
    newdat$level <- level
    newdat$response <- response
    newdat
  })
  
  bind_rows(out)
}

# ---- Generate predictions ----
preds_elev_iWUE  <- make_preds_with_species(gam_iWUE_Spe, raw.env.data, "Elevation",      "iWUE", target_species)
preds_T2_iWUE    <- make_preds_with_species(gam_iWUE_Spe, raw.env.data, "mean_T2",        "iWUE", target_species)
preds_moist_iWUE <- make_preds_with_species(gam_iWUE_Spe, raw.env.data, "mean_moist_pct", "iWUE", target_species)

# ---- Plotting function ----
purple_cols <- c("Low"="#d3bdf0", "Medium"="#9b59b6", "High"="#5e3370")

plot_driver_dummy <- function(preds, driver_label, data, driver, target_species) {
  response <- unique(preds$response)
  qs <- quantile(data[[driver]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
  names(qs) <- c("Low", "Medium", "High")
  
  data$level <- cut(data[[driver]],
                    breaks = c(-Inf, qs["Low"], qs["Medium"], qs["High"], Inf),
                    labels = c("Low", "Medium", "High", NA),
                    include.lowest = TRUE, right = TRUE)
  
  ggplot(preds, aes(x = Tleaf, y = fit, color = level)) +
    geom_point(data = data %>% filter(Species == target_species),
               aes(x = Tleaf, y = .data[[response]], color = level),
               inherit.aes = FALSE, alpha = 0.01, size = 0.5) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se, fill = level),
                alpha = 0.2, colour = NA) +
    scale_color_manual(values = purple_cols, breaks = c("Low", "Medium", "High")) +
    scale_fill_manual(values = purple_cols, breaks = c("Low", "Medium", "High")) +
    theme_classic() +
    coord_cartesian(ylim = c(0, 250)) +
    labs(x = "Tleaf [°C]", y = "iWUE [µmol mol⁻¹]",
         color = driver_label, fill = driver_label,
         title = paste("Predicted iWUE for", target_species))
}

# ---- Final three-panel plot ----
p1_iWUE <- plot_driver_dummy(preds_elev_iWUE,  "Elevation",   raw.env.data, "Elevation",      target_species)
p2_iWUE <- plot_driver_dummy(preds_T2_iWUE,    "Temperature", raw.env.data, "mean_T2",        target_species)
p3_iWUE <- plot_driver_dummy(preds_moist_iWUE, "Moisture",    raw.env.data, "mean_moist_pct", target_species)

ggarrange(p1_iWUE, p2_iWUE, p3_iWUE, nrow = 1, labels = c("A", "B", "C"))




#


























































# ---- Fit GAM with species dummy variables ----
species_levels <- levels(raw.env.data$Species)
ref_species <- species_levels[1]
other_species <- species_levels[-1]

raw2 <- raw.env.data %>%
  mutate(across(Species, as.character))

# Create dummy variables
for(sp in other_species) {
  nm <- make.names(sp)
  raw2[[nm]] <- as.integer(raw2$Species == sp)
}

# Build formula
intercept_terms <- paste0("`", make.names(other_species), "`", collapse = " + ")
smooth_terms <- paste0("s(Tleaf, by = `", make.names(other_species), "`, k = 4)", collapse = " + ")

fml <- as.formula(paste(
  "iWUE ~ s(Tleaf, k=4) +",
  "s(Elevation, k=4) + s(mean_T2, k=4) + s(mean_moist_pct, k=4) +",
  intercept_terms, "+",
  smooth_terms, "+",
  "ti(Tleaf, Elevation, k=4) + ti(Tleaf, mean_T2, k=4)"
))

gam_dummy <- gam(fml, data = raw2, method = "REML")


ref_species <- names(sort(table(raw.env.data$Species), decreasing = TRUE))[1]

make_preds_with_species <- function(model, data, driver, response, ref_species) {
  qs <- quantile(data[[driver]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
  names(qs) <- c("Low", "Medium", "High")
  
  out <- lapply(names(qs), function(level) {
    newdat <- expand.grid(
      Tleaf = seq(min(data$Tleaf, na.rm = TRUE),
                  max(data$Tleaf, na.rm = TRUE), length.out = 200),
      Elevation = median(data$Elevation, na.rm = TRUE),
      mean_T2   = median(data$mean_T2, na.rm = TRUE),
      mean_moist_pct = median(data$mean_moist_pct, na.rm = TRUE),
      Species = ref_species
    )
    newdat[[driver]] <- qs[level]
    
    p <- predict(model, newdata = newdat, se.fit = TRUE)
    newdat$fit <- p$fit
    newdat$se <- p$se.fit
    newdat$driver <- driver
    newdat$level <- level
    newdat$response <- response
    newdat
  })
  
  bind_rows(out)
}

preds_elev_iWUE  <- make_preds_with_species(gam_iWUE_Spe, raw.env.data, "Elevation",      "iWUE", ref_species)
preds_T2_iWUE    <- make_preds_with_species(gam_iWUE_Spe, raw.env.data, "mean_T2",        "iWUE", ref_species)
preds_moist_iWUE <- make_preds_with_species(gam_iWUE_Spe, raw.env.data, "mean_moist_pct", "iWUE", ref_species)
# ---- Plotting function ----
purple_cols <- c("Low"="#d3bdf0", "Medium"="#9b59b6", "High"="#5e3370")

plot_driver_dummy <- function(preds, driver_label, data, driver) {
  response <- unique(preds$response)
  qs <- quantile(data[[driver]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
  names(qs) <- c("Low", "Medium", "High")
  
  data$level <- cut(data[[driver]],
                    breaks = c(-Inf, qs["Low"], qs["Medium"], qs["High"], Inf),
                    labels = c("Low", "Medium", "High", NA),
                    include.lowest = TRUE, right = TRUE)
  
  ggplot(preds, aes(x = Tleaf, y = fit, color = level)) +
    geom_point(data = data %>% filter(Species == target_species),
               aes(x = Tleaf, y = .data[[response]], color = level),
               inherit.aes = FALSE, alpha = 0.01, size = 0.5) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se, fill = level),
                alpha = 0.2, colour = NA) +
    scale_color_manual(values = purple_cols, breaks = c("Low", "Medium", "High")) +
    scale_fill_manual(values = purple_cols, breaks = c("Low", "Medium", "High")) +
    theme_classic() +
    labs(x = "Tleaf [°C]", y = "iWUE [µmol mol⁻¹]",
         color = driver_label, fill = driver_label,
         title = paste("Predicted iWUE for", target_species))
}


p1_iWUE <- plot_driver_dummy(preds_elev_iWUE,  "Elevation",   raw.env.data, "Elevation")
p2_iWUE <- plot_driver_dummy(preds_T2_iWUE,    "Temperature", raw.env.data, "mean_T2")
p3_iWUE <- plot_driver_dummy(preds_moist_iWUE, "Moisture",    raw.env.data, "mean_moist_pct")

ggarrange(p1_iWUE, p2_iWUE, p3_iWUE, nrow = 1, labels = c("A", "B", "C"))



















#Producing Whack lines

# ---- Prediction function for one species ----
make_preds_species <- function(model, data, driver, response,
                               species_levels, other_species, target_species) {
  qs <- quantile(data[[driver]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
  names(qs) <- c("Low", "Medium", "High")
  
  out <- lapply(names(qs), function(level) {
    newdat <- expand.grid(
      Tleaf = seq(min(data$Tleaf, na.rm = TRUE),
                  max(data$Tleaf, na.rm = TRUE), length.out = 200),
      Elevation      = median(data$Elevation, na.rm = TRUE),
      mean_T2        = median(data$mean_T2, na.rm = TRUE),
      mean_moist_pct = median(data$mean_moist_pct, na.rm = TRUE)
    )
    
    # Set species dummies
    for(sp in other_species) {
      newdat[[make.names(sp)]] <- as.integer(sp == target_species)
    }
    
    newdat[[driver]] <- qs[level]
    
    p <- predict(model, newdata = newdat, se.fit = TRUE)
    newdat$fit <- p$fit
    newdat$se <- p$se.fit
    newdat$driver <- driver
    newdat$level <- level
    newdat$response <- response
    newdat$Species <- target_species
    newdat
  })
  
  bind_rows(out)
}

# ---- Generate predictions for one species ----
target_species <- "Senecio_glaberrimus"  # Change this to any species you want

preds_elev_iWUE  <- make_preds_species(gam_dummy, raw2, "Elevation",      "iWUE", species_levels, other_species, target_species)
preds_T2_iWUE    <- make_preds_species(gam_dummy, raw2, "mean_T2",        "iWUE", species_levels, other_species, target_species)
preds_moist_iWUE <- make_preds_species(gam_dummy, raw2, "mean_moist_pct", "iWUE", species_levels, other_species, target_species)


# ---- Final three-panel figure ----
p1_iWUE <- plot_driver_dummy(preds_elev_iWUE,  "Elevation",   raw2, "Elevation")
p2_iWUE <- plot_driver_dummy(preds_T2_iWUE,    "Temperature", raw2, "mean_T2")
p3_iWUE <- plot_driver_dummy(preds_moist_iWUE, "Moisture",    raw2, "mean_moist_pct")

ggarrange(p1_iWUE, p2_iWUE, p3_iWUE, nrow = 1, labels = c("A", "B", "C"))

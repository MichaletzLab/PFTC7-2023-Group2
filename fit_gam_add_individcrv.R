# Function to fit GAM models
fit_gam <- function(data, response, predictors) {
  gam_formula <- as.formula(paste(response, "~", paste(sapply(predictors, function(x) paste("s(", x, ", bs='ts', k=5)")), collapse = "+")))
  fit <- gam(gam_formula, data = data, select = TRUE, method = "REML")
  return(fit)
}

# Fit both GAM models (tleaf only and tleaf + cond) for each group or condition
curve_ids <- unique(dat$curveid)
gam_models <- lapply(curve_ids, function(curve) {
  data_subset <- dat[curveid == curve]
  model_tleaf <- fit_gam(data_subset, "photo", "tleaf")
  model_tleaf_cond <- fit_gam(data_subset, "photo", c("tleaf", "cond"))
  list(tleaf = model_tleaf, tleaf_cond = model_tleaf_cond)
})

# Extract smooth estimates and SE for tleaf for each group and model
smooth_estimates <- lapply(gam_models, function(models) {
  tleaf_estimates <- gratia::smooth_estimates(models$tleaf, unconditional = TRUE, smooth = 's(tleaf)')
  tleaf_cond_estimates <- gratia::smooth_estimates(models$tleaf_cond, unconditional = TRUE, smooth = 's(tleaf)')
  list(tleaf = tleaf_estimates, tleaf_cond = tleaf_cond_estimates)
})

# Ensure consistent column names and combine smooth estimates into a single data table with consistent column names
smooth_data <- rbindlist(lapply(seq_along(curve_ids), function(i) {
  curveid <- curve_ids[i]
  tleaf_data <- data.table(curveid = curveid, model = "tleaf", smooth_estimates[[i]]$tleaf)
  tleaf_cond_data <- data.table(curveid = curveid, model = "tleaf_cond", smooth_estimates[[i]]$tleaf_cond)
  
  # Ensure consistent column names
  tleaf_data <- tleaf_data[, .(curveid, model, tleaf, .estimate, .se)]
  tleaf_cond_data <- tleaf_cond_data[, .(curveid, model, tleaf, .estimate, .se)]
  
  rbind(tleaf_data, tleaf_cond_data, fill = TRUE)
}), fill = TRUE)

# Convert curveid and model to factor
smooth_data[, curveid := factor(curveid, levels = curve_ids)]
smooth_data[, model := factor(model, levels = c("tleaf", "tleaf_cond"))]

# Function to find the optimal temperature (tleaf at maxA)
find_optimal_tleaf <- function(smooth_data) {
  optimal_row <- smooth_data[which.max(.estimate)]
  return(optimal_row$tleaf)
}

# Function to calculate breadth
get_breadth <- function(photo, tleaf, threshold = 0.95) {
  # Find the indices where photo drops below the threshold
  idx_left <- max(which(photo >= threshold * max(photo)), default = 1)
  idx_right <- min(which(photo >= threshold * max(photo)), default = length(photo))
  
  # Calculate breadth as absolute difference
  breadth <- abs(tleaf[idx_right] - tleaf[idx_left])
  
  return(list(idx_left = idx_left, idx_right = idx_right, breadth = breadth))
}

# Data frame to store optimal temperatures and breadth
optimal_temps <- data.table(curveid = character(),
                            tleaf_opt = numeric(),
                            tleaf_cond_opt = numeric(),
                            weib_Topt = numeric(),
                            weib_breadth = numeric(),
                            tleaf_breadth = numeric(),
                            tleaf_cond_breadth = numeric())

# Loop through each curveid
for (i in seq_along(curve_ids)) {
  curve_data <- dat[curveid == curve_ids[i]]
  smooth_data_curve <- smooth_data[curveid == curve_ids[i]]
  
  # Calculate breadths and get indices using get_breadth function
  tleaf_breadth_data <- get_breadth(smooth_data_curve[model == "tleaf"]$.estimate, smooth_data_curve[model == "tleaf"]$tleaf)
  tleaf_cond_breadth_data <- get_breadth(smooth_data_curve[model == "tleaf_cond"]$.estimate, smooth_data_curve[model == "tleaf_cond"]$tleaf)
  
  # Extract breadths
  tleaf_breadth <- tleaf_breadth_data$breadth
  tleaf_cond_breadth <- tleaf_cond_breadth_data$breadth
  
  # Append optimal temperatures and breadths to the data frame
  optimal_temps <- rbind(optimal_temps, data.table(
    curveid = curve_ids[i],
    tleaf_opt = optimal_tleaf_tleaf,
    tleaf_cond_opt = optimal_tleaf_tleaf_cond,
    weib_Topt = weib_Topt_value,
    weib_breadth = weib_breadth_value,
    tleaf_breadth = tleaf_breadth,
    tleaf_cond_breadth = tleaf_cond_breadth
  ))
}


# Plotting each individual curve with both GAMs and optimal temperature lines
plots <- lapply(seq_along(curve_ids), function(i) {
  curve_data <- dat[curveid == curve_ids[i]]
  smooth_data_curve <- smooth_data[curveid == curve_ids[i]]
  
  # Extract unique weib.topt and weib.breadth values for the current curveid
  weib_Topt_value <- unique(curve_data$weib.topt)
  weib_breadth_value <- unique(curve_data$weib.breadth)
  
  optimal_tleaf_tleaf <- find_optimal_tleaf(smooth_data_curve[model == "tleaf"])
  optimal_tleaf_tleaf_cond <- find_optimal_tleaf(smooth_data_curve[model == "tleaf_cond"])
  
  # Calculate breadths
  tleaf_breadth <- get_breadth(smooth_data_curve[model == "tleaf"]$.estimate, smooth_data_curve[model == "tleaf"]$tleaf, max(smooth_data_curve[model == "tleaf"]$.estimate))
  tleaf_cond_breadth <- get_breadth(smooth_data_curve[model == "tleaf_cond"]$.estimate, smooth_data_curve[model == "tleaf_cond"]$tleaf, max(smooth_data_curve[model == "tleaf_cond"]$.estimate))
  
  # Append optimal temperatures and breadths to the data frame
  optimal_temps <<- rbind(optimal_temps, data.table(curveid = curve_ids[i], tleaf_opt = optimal_tleaf_tleaf, tleaf_cond_opt = optimal_tleaf_tleaf_cond, weib_Topt = weib_Topt_value, weib_breadth = weib_breadth_value, tleaf_breadth = tleaf_breadth, tleaf_cond_breadth = tleaf_cond_breadth))
  
  # Get the threshold * peak value for each model
  threshold_peak_value_tleaf <- 0.95 * max(smooth_data_curve[model == "tleaf"]$.estimate)
  threshold_peak_value_tleaf_cond <- 0.95 * max(smooth_data_curve[model == "tleaf_cond"]$.estimate)
  
  # Find the indices where the curve crosses the threshold peak value
  idx_cross_tleaf <- which(smooth_data_curve[model == "tleaf"]$.estimate >= threshold_peak_value_tleaf)
  idx_cross_tleaf_cond <- which(smooth_data_curve[model == "tleaf_cond"]$.estimate >= threshold_peak_value_tleaf_cond)
  
  # Calculate the x-values for the horizontal lines
  x_start_tleaf <- min(smooth_data_curve[model == "tleaf"]$tleaf[idx_cross_tleaf])
  x_end_tleaf <- max(smooth_data_curve[model == "tleaf"]$tleaf[idx_cross_tleaf])
  x_start_tleaf_cond <- min(smooth_data_curve[model == "tleaf_cond"]$tleaf[idx_cross_tleaf_cond])
  x_end_tleaf_cond <- max(smooth_data_curve[model == "tleaf_cond"]$tleaf[idx_cross_tleaf_cond])
  
  ggplot() +
    geom_line(data = smooth_data_curve[model == "tleaf"], 
              aes(x = tleaf, y = .estimate, color = "tleaf"), size = 1.2) +
    geom_line(data = smooth_data_curve[model == "tleaf_cond"], 
              aes(x = tleaf, y = .estimate, color = "tleaf + cond"), size = 1.2) +
    geom_ribbon(data = smooth_data_curve[model == "tleaf"], 
                aes(x = tleaf, ymin = .estimate - 1.96 * .se, ymax = .estimate + 1.96 * .se, fill = "tleaf"), alpha = 0.2) +
    geom_ribbon(data = smooth_data_curve[model == "tleaf_cond"], 
                aes(x = tleaf, ymin = .estimate - 1.96 * .se, ymax = .estimate + 1.96 * .se, fill = "tleaf + cond"), alpha = 0.2) +
    geom_vline(xintercept = optimal_tleaf_tleaf, color = "blue", linetype = "dashed", size = 1) +
    geom_vline(xintercept = optimal_tleaf_tleaf_cond, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = weib_Topt_value, color = "grey", linetype = "dashed", size = 1) +
    geom_segment(aes(x = x_start_tleaf, y = threshold_peak_value_tleaf,
                     xend = x_end_tleaf, yend = threshold_peak_value_tleaf),
                 color = "blue", size = 1.2) +
    geom_segment(aes(x = x_start_tleaf_cond, y = threshold_peak_value_tleaf_cond,
                     xend = x_end_tleaf_cond, yend = threshold_peak_value_tleaf_cond),
                 color = "red", size = 1.2) +
    labs(title = paste("Curve:", curve_ids[i]), x = "Leaf Temperature (°C)", y = "Partial Effect on Photosynthesis (µmol m^-2 s^-1)") +
    scale_color_manual(name = "Model", 
                       values = c("tleaf" = "blue", "tleaf + cond" = "red"),
                       labels = c("tleaf" = "tleaf", "tleaf + cond" = "tleaf + cond")) +
    scale_fill_manual(name = "Model", 
                      values = c("tleaf" = "blue", "tleaf + cond" = "red"),
                      labels = c("tleaf" = "tleaf", "tleaf + cond" = "tleaf + cond")) +
    theme_minimal()
})

# Save the plots to a PDF with one plot per page
pdf("gam_plots_with_optimal_temp.pdf")
for (plot in plots) {
  print(plot)
}
dev.off()
fwrite(optimal_temps, "gam.fits.add.topts.csv", row.names = FALSE)

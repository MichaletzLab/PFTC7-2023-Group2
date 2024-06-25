# Plotting each individual curve with both GAMs and optimal temperature lines
plots <- lapply(seq_along(curve_ids), function(i) {
  curve_data <- dat[curveid == curve_ids[i]]
  smooth_data_curve <- smooth_data[curveid == curve_ids[i]]
  
  # Extract unique weib.topt and weib.breadth values for the current curveid
  weib_Topt_value <- unique(curve_data$weib.topt)
  weib_breadth_value <- unique(curve_data$weib.breadth)
  
  # Calculate optimal temperatures
  optimal_tleaf_tleaf <- find_optimal_tleaf(smooth_data_curve[model == "tleaf"])
  optimal_tleaf_tleaf_cond <- find_optimal_tleaf(smooth_data_curve[model == "tleaf_cond"])
  
  # Append optimal temperatures to optimal_temps data table
  optimal_temps <- rbind(optimal_temps, data.table(
    curveid = curve_ids[i],
    tleaf_opt = optimal_tleaf_tleaf,
    tleaf_cond_opt = optimal_tleaf_tleaf_cond,
    weib_Topt = weib_Topt_value,
    weib_breadth = weib_breadth_value,
    tleaf_breadth = NA,  # Update this with actual calculations
    tleaf_cond_breadth = NA  # Update this with actual calculations
  ))
  
  # Calculate breadths
  tleaf_breadth <- get_breadth(smooth_data_curve[model == "tleaf"]$.estimate, smooth_data_curve[model == "tleaf"]$tleaf, max(smooth_data_curve[model == "tleaf"]$.estimate))
  tleaf_cond_breadth <- get_breadth(smooth_data_curve[model == "tleaf_cond"]$.estimate, smooth_data_curve[model == "tleaf_cond"]$tleaf, max(smooth_data_curve[model == "tleaf_cond"]$.estimate))
  
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

# Write optimal_temps to CSV
fwrite(optimal_temps, "gam.fits.add.topts.csv", row.names = FALSE)


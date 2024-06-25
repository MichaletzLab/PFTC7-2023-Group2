# Function to calculate breadth from predictions
get_breadth <- function(predicted_photo, tleaf, threshold = 0.95) {
  max_photo <- max(predicted_photo, na.rm = TRUE)
  threshold_value <- threshold * max_photo
  
  idx_left <- which(predicted_photo >= threshold_value)[1]
  idx_right <- which(predicted_photo >= threshold_value)[length(which(predicted_photo >= threshold_value))]
  
  # Calculate breadth as absolute difference
  breadth <- abs(tleaf[idx_right] - tleaf[idx_left])
  
  return(list(idx_left = idx_left, idx_right = idx_right, breadth = breadth))
}

# Function to fit GLMs and plot for each curveid
fit_glm_and_plot <- function(data, 
                             curveid_col = "curveid", 
                             x_col = "tleaf", 
                             y_col = "photo", 
                             cond_col = "cond", 
                             topt_col = "weib.topt", 
                             breadth_col = "weib.breadth") {
  
  # Initialize results dataframe
  results <- data.frame(curveid = numeric(), 
                        topt_temp = numeric(), 
                        topt_temp_cond = numeric(), 
                        topt_weib = numeric(),
                        breadth_temp = numeric(), 
                        breadth_temp_cond = numeric(), 
                        breadth_weib = numeric())
  
  # Create a PDF file to save the plots
  pdf("glm_fits_curves.pdf", width = 8, height = 8)
  
  # Get unique curve IDs
  unique_curve_ids <- unique(data[[curveid_col]])
  
  # Loop over each curve ID
  for (curve_id in unique_curve_ids) {
    curve_data <- subset(data, data[[curveid_col]] == curve_id)
    print(paste("Fitting curve ID:", curve_id))
    
    # Fit GLM with temperature and conductance
    glm_temp_cond <- tryCatch(
      glm(photo ~ tleaf + I(tleaf^2) + cond, data = curve_data), 
      error = function(e) {
        print(paste("BadFit: glm_temp_cond for Curve ID", curve_id))
        print(e)
        return(NULL)
      }
    )
    
    # Fit GLM with only temperature
    glm_temp <- tryCatch(
      glm(photo ~ tleaf + I(tleaf^2), data = curve_data), 
      error = function(e) {
        print(paste("BadFit: glm_temp for Curve ID", curve_id))
        print(e)
        return(NULL)
      }
    )
    
    # Check if both models were successfully fitted
    if (!is.null(glm_temp_cond) && !is.null(glm_temp)) {
      # Predictions and Topt calculation for temperature + conductance model
      curve_data$pred_temp_cond <- predict(glm_temp_cond, newdata = curve_data)
      topt_temp_cond <- curve_data$tleaf[which.max(curve_data$pred_temp_cond)]
      breadth_temp_cond <- get_breadth(curve_data$pred_temp_cond, curve_data$tleaf, threshold = 0.95)
      
      # Predictions and Topt calculation for temperature model
      curve_data$pred_temp <- predict(glm_temp, newdata = curve_data)
      topt_temp <- curve_data$tleaf[which.max(curve_data$pred_temp)]
      breadth_temp <- get_breadth(curve_data$pred_temp, curve_data$tleaf, threshold = 0.95)
      
      # Calculate Weibull Topt and breadth
      topt_weib <- mean(curve_data[[topt_col]], na.rm = TRUE)
      breadth_weib <- mean(curve_data[[breadth_col]], na.rm = TRUE)
      
      # Store results
      results <- rbind(results, data.frame(curveid = curve_id, 
                                           topt_temp = topt_temp, 
                                           topt_temp_cond = topt_temp_cond, 
                                           topt_weib = topt_weib,
                                           breadth_temp = breadth_temp$breadth,
                                           breadth_temp_cond = breadth_temp_cond$breadth,
                                           breadth_weib = breadth_weib))
      
      # Create the plot
      plot <- ggplot(curve_data, aes_string(x = x_col, y = y_col)) +
        geom_point() +
        geom_line(aes(y = pred_temp), color = 'red') +
        geom_line(aes(y = pred_temp_cond), color = 'blue') +
        geom_vline(xintercept = topt_temp, color = 'red', linetype = "dashed") +
        geom_vline(xintercept = topt_temp_cond, color = 'blue', linetype = "dashed") +
        geom_vline(xintercept = topt_weib, color = 'grey', linetype = "dashed") +
        geom_segment(aes(x = curve_data[[x_col]][breadth_temp$idx_left], 
                         y = max(curve_data$pred_temp) * 0.95,
                         xend = curve_data[[x_col]][breadth_temp$idx_right], 
                         yend = max(curve_data$pred_temp) * 0.95),
                     color = "red", size = 1.2) +
        geom_segment(aes(x = curve_data[[x_col]][breadth_temp_cond$idx_left], 
                         y = max(curve_data$pred_temp_cond) * 0.95,
                         xend = curve_data[[x_col]][breadth_temp_cond$idx_right], 
                         yend = max(curve_data$pred_temp_cond) * 0.95),
                     color = "blue", size = 1.2) +
        theme_classic() +
        ggtitle(paste("Curve ID:", curve_id))
      
      # Print plot to PDF
      print(plot)
    }
  }
  
  # Close the PDF device
  dev.off()
  
  return(results)
}

# Call the function with your data
results <- fit_glm_and_plot(dat)
write.csv(results, "glm.fits.topts.csv")

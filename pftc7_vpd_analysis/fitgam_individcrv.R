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

# Function to fit GAM models and plot, calculate AIC, and concurvity
fit_gams_and_plot <- function(data, curveid_col = "curveid", x_col = "tleaf", y_col = "photo", cond_col = "cond", topt_col = "weib.topt", breadth_col = "weib.breadth") {
  
  models <- c("gam_temp_cond", "gam_temp")
  results <- data.frame(curveid = numeric(), topt_temp = numeric(), topt_temp_cond = numeric(), topt_weib = numeric(), breadth_temp = numeric(), breadth_temp_cond = numeric(), breadth_weib = numeric(), AIC_tleaf = numeric(), AIC_tleaf_cond = numeric(), concurvity_tleaf = numeric(), concurvity_tleaf_cond = numeric())
  
  # Create a PDF file to save the plots
  pdf("gam_fits.pdf", width = 8, height = 8)  # Adjust width and height to make it more square
  
  for (curve_id in unique(data[[curveid_col]])) {
    curve_data <- subset(data, data[[curveid_col]] == curve_id) # Subset data for the current curveID
    print(paste("Curve ID:", curve_id, ", Number of observations:", nrow(curve_data)))  # Debug
    
    # Initialize prediction columns
    curve_data$pred_temp <- NA
    curve_data$pred_temp_cond <- NA
    
    topt_values <- list(topt_temp = NA, topt_temp_cond = NA)
    breadth_values <- list(breadth_temp = NA, breadth_temp_cond = NA)
    AIC_values <- list(AIC_tleaf = NA, AIC_tleaf_cond = NA)
    concurvity_values <- list(concurvity_tleaf = NA, concurvity_tleaf_cond = NA)
    
    for (model in models) {
      fit <- NULL
      if (model == "gam_temp_cond") {
        # Fit GAM with temperature and conductance
        fit <- tryCatch(
          gam(photo ~ s(cond, k = 5, bs = 'ts') + s(tleaf, k = 5, bs = 'ts'), 
              data = curve_data,
              select = TRUE,
              method = "REML"), 
          error = function(e) {
            print("BadFit: gam_temp_cond")
            return(NULL)
          }
        )
        if (!is.null(fit)) {
          curve_data$pred_temp_cond <- predict(fit, newdata = curve_data)
          topt_values$topt_temp_cond <- gratia::smooth_estimates(fit, unconditional = TRUE, smooth = 's(tleaf)') %>% 
            filter(.estimate == max(.estimate)) %>% 
            pull(tleaf)
          breadth_values_temp_cond <- get_breadth(curve_data$pred_temp_cond, curve_data$tleaf, threshold=0.95)
          breadth_values$breadth_temp_cond <- breadth_values_temp_cond$breadth
          AIC_values$AIC_tleaf_cond <- AIC(fit)
          concurvity_values$concurvity_tleaf_cond <- concurvity(fit)
          
          # Calculate start and end points for the horizontal line
          x_start_tleaf_cond <- curve_data$tleaf[breadth_values_temp_cond$idx_left]
          x_end_tleaf_cond <- curve_data$tleaf[breadth_values_temp_cond$idx_right]
          y_value_tleaf_cond <- max(curve_data$pred_temp_cond) * 0.95
        }
      } else if (model == "gam_temp") {
        # Fit GAM with only temperature
        fit <- tryCatch(
          gam(photo ~ s(tleaf, bs='ts', k=5), 
              data = curve_data,
              select = TRUE,
              method = "REML"), 
          error = function(e) {
            print("BadFit: gam_temp")
            return(NULL)
          }
        )
        if (!is.null(fit)) {
          curve_data$pred_temp <- predict(fit, newdata = curve_data)
          topt_values$topt_temp <- gratia::smooth_estimates(fit, unconditional = TRUE, smooth = 's(tleaf)') %>% 
            filter(.estimate == max(.estimate)) %>% 
            pull(tleaf)
          breadth_values_temp <- get_breadth(curve_data$pred_temp, curve_data$tleaf, threshold=0.95)
          breadth_values$breadth_temp <- breadth_values_temp$breadth
          AIC_values$AIC_tleaf <- AIC(fit)
          concurvity_values$concurvity_tleaf <- concurvity(fit)
          
          # Calculate start and end points for the horizontal line
          x_start_tleaf <- curve_data$tleaf[breadth_values_temp$idx_left]
          x_end_tleaf <- curve_data$tleaf[breadth_values_temp$idx_right]
          y_value_tleaf <- max(curve_data$pred_temp) * 0.95
        }
      }
    }
    
    topt_weib <- mean(curve_data[[topt_col]], na.rm = TRUE)
    topt_values$topt_weib <- topt_weib
    breadth_weib <- mean(curve_data[[breadth_col]], na.rm = TRUE)
    breadth_values$breadth_weib <- breadth_weib
    
    # Save the Topt, breadth, AIC, and concurvity values in the results dataframe
    results <- rbind(results, data.frame(curveid = curve_id, 
                                         topt_temp = topt_values$topt_temp, 
                                         topt_temp_cond = topt_values$topt_temp_cond, 
                                         topt_weib = topt_weib,
                                         breadth_temp = breadth_values$breadth_temp,
                                         breadth_temp_cond = breadth_values$breadth_temp_cond,
                                         breadth_weib = breadth_weib,
                                         AIC_tleaf = AIC_values$AIC_tleaf,
                                         AIC_tleaf_cond = AIC_values$AIC_tleaf_cond,
                                         concurvity_tleaf = concurvity_values$concurvity_tleaf,
                                         concurvity_tleaf_cond = concurvity_values$concurvity_tleaf_cond))
    
    # Initialize plot
    plot <- ggplot(curve_data, aes(x = tleaf, y = photo)) +
      geom_point() +
      theme_classic() +
      ggtitle(paste("Curve ID:", curve_id))
    
    # Add prediction lines
    if (!is.na(curve_data$pred_temp[1])) {
      plot <- plot + geom_line(aes(y = pred_temp), color = 'red')
    }
    if (!is.na(curve_data$pred_temp_cond[1])) {
      plot <- plot + geom_line(aes(y = pred_temp_cond), color = 'blue')
    }
    
    # Add vertical lines for Topt values
    if (!is.na(topt_values$topt_temp)) {
      plot <- plot + geom_vline(xintercept = topt_values$topt_temp, color = 'red', linetype = "dashed")
    }
    if (!is.na(topt_values$topt_temp_cond)) {
      plot <- plot + geom_vline(xintercept = topt_values$topt_temp_cond, color = 'blue', linetype = "dashed")
    }
    if (!is.na(topt_weib)) {
      plot <- plot + geom_vline(xintercept = topt_weib, color = 'grey', linetype = "dashed")
    }
    
    # Add horizontal lines at the correct y-value
    if (!is.na(breadth_values$breadth_temp_cond)) {
      plot <- plot + geom_segment(aes(x = x_start_tleaf_cond, y = y_value_tleaf_cond,
                                      xend = x_end_tleaf_cond, yend = y_value_tleaf_cond),
                                  color = "blue", size = 1.2)
    }
    if (!is.na(breadth_values$breadth_temp)) {
      plot <- plot + geom_segment(aes(x = x_start_tleaf, y = y_value_tleaf,
                                      xend = x_end_tleaf, yend = y_value_tleaf),
                                  color = "red", size = 1.2)
    }
    
    # Add annotation for legend in top right corner dynamically
    plot <- plot +
      annotation_custom(
        grob = textGrob("Temperature GAM", gp = gpar(col = "red", fontsize = 8)),
        xmin = max(curve_data$tleaf) * 0.8, ymin = max(curve_data$photo) * 0.95
      ) +
      annotation_custom(
        grob = textGrob("Temperature + Conductance GAM", gp = gpar(col = "blue", fontsize = 8)),
        xmin = max(curve_data$tleaf) * 0.8, ymin = max(curve_data$photo) * 0.9
      )
    # Add legend
    plot <- plot +
      labs(x = "Temperature", y = "Photo") +
      theme(legend.position = c(0.8, 0.2))
    
    # Save plot to PDF
    print(plot)
    
  }
  
  # Close the PDF device
  dev.off()
  
  # Return results dataframe
  return(results)
}

# Example usage:
results <- fit_gams_and_plot(dat)
filtered_results <- results %>%
  rownames_to_column(var = "row_label")%>%
  filter(str_detect(row_label, "worst"))

write.csv(filtered_results, "gam.fits.raw.topts.csv")

    
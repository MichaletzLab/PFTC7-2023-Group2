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
fit_gams_and_plot <- function(data, curveid_col = "curveid", x_col = "tleaf", y_col = "photo", topt_col = "weib.topt", breadth_col = "weib.breadth") {
  
  models <- c("gam_tleaf_cond", "gam_tleaf", "gam_tleaf_vpd", "gam_tleaf_vpd_cond")
  results <- data.frame(curveid = numeric(), topt_tleaf = numeric(), topt_tleaf_cond = numeric(), topt_tleaf_vpdl = numeric(), topt_tleaf_vpdl_cond = numeric(), topt_weib = numeric(), breadth_tleaf = numeric(), breadth_tleaf_cond = numeric(), breadth_tleaf_vpdl = numeric(), breadth_tleaf_vpdl_cond = numeric(), breadth_weib = numeric(), AIC_tleaf = numeric(), AIC_tleaf_cond = numeric(), AIC_tleaf_vpdl = numeric(), AIC_tleaf_vpdl_cond = numeric(), concurvity_tleaf = numeric(), concurvity_tleaf_cond = numeric(), concurvity_tleaf_vpdl = numeric(), concurvity_tleaf_vpdl_cond = numeric())
  
  # Create a PDF file to save the plots
  pdf("gam_fits.pdf", width = 8, height = 8)  # Adjust width and height to make it more square
  
  for (curve_id in unique(data[[curveid_col]])) {
    curve_data <- subset(data, data[[curveid_col]] == curve_id) # Subset data for the current curveID
    print(paste("Curve ID:", curve_id, ", Number of observations:", nrow(curve_data)))  # Debug
    
    # Initialize prediction columns
    curve_data$pred_tleaf <- NA
    curve_data$pred_tleaf_cond <- NA
    curve_data$pred_tleaf_vpdl <- NA
    curve_data$pred_tleaf_vpdl_cond <- NA
    
    topt_values <- list(topt_tleaf = NA, topt_tleaf_cond = NA, topt_tleaf_vpdl = NA, topt_tleaf_vpdl_cond = NA)
    breadth_values <- list(breadth_tleaf = NA, breadth_tleaf_cond = NA, breadth_tleaf_vpdl = NA, breadth_tleaf_vpdl_cond = NA)
    AIC_values <- list(AIC_tleaf = NA, AIC_tleaf_cond = NA, AIC_tleaf_vpdl = NA, AIC_tleaf_vpdl_cond = NA)
    concurvity_values <- list(concurvity_tleaf = NA, concurvity_tleaf_cond = NA, concurvity_tleaf_vpdl = NA, concurvity_tleaf_vpdl_cond = NA)
    
    for (model in models) {
      fit <- NULL
      if (model == "gam_tleaf_cond") {
        # Fit GAM with tleaferature and conductance
        fit <- tryCatch(
          gam(photo ~ s(cond, k = 5, bs = 'ts') + s(tleaf, k = 5, bs = 'ts'), 
              data = curve_data,
              select = TRUE,
              method = "REML"), 
          error = function(e) {
            print("BadFit: gam_tleaf_cond")
            return(NULL)
          }
        )
        if (!is.null(fit)) {
          curve_data$pred_tleaf_cond <- predict(fit, newdata = curve_data)
          topt_values$topt_tleaf_cond <- gratia::smooth_estimates(fit, unconditional = TRUE, smooth = 's(tleaf)') %>% 
            filter(.estimate == max(.estimate)) %>% 
            pull(tleaf)
          breadth_values_tleaf_cond <- get_breadth(curve_data$pred_tleaf_cond, curve_data$tleaf, threshold=0.95)
          breadth_values$breadth_tleaf_cond <- breadth_values_tleaf_cond$breadth
          AIC_values$AIC_tleaf_cond <- AIC(fit)
          concurvity_values$concurvity_tleaf_cond <- concurvity(fit)
          
          # Calculate start and end points for the horizontal line
          x_start_tleaf_cond <- curve_data$tleaf[breadth_values_tleaf_cond$idx_left]
          x_end_tleaf_cond <- curve_data$tleaf[breadth_values_tleaf_cond$idx_right]
          y_value_tleaf_cond <- max(curve_data$pred_tleaf_cond) * 0.95
        }
      } else if (model == "gam_tleaf") {
        # Fit GAM with only tleaferature
        fit <- tryCatch(
          gam(photo ~ s(tleaf, bs='ts', k=5), 
              data = curve_data,
              select = TRUE,
              method = "REML"), 
          error = function(e) {
            print("BadFit: gam_tleaf")
            return(NULL)
          }
        )
        if (!is.null(fit)) {
          curve_data$pred_tleaf <- predict(fit, newdata = curve_data)
          topt_values$topt_tleaf <- gratia::smooth_estimates(fit, unconditional = TRUE, smooth = 's(tleaf)') %>% 
            filter(.estimate == max(.estimate)) %>% 
            pull(tleaf)
          breadth_values_tleaf <- get_breadth(curve_data$pred_tleaf, curve_data$tleaf, threshold=0.95)
          breadth_values$breadth_tleaf <- breadth_values_tleaf$breadth
          AIC_values$AIC_tleaf <- AIC(fit)
          concurvity_values$concurvity_tleaf <- concurvity(fit)
          
          # Calculate start and end points for the horizontal line
          x_start_tleaf <- curve_data$tleaf[breadth_values_tleaf$idx_left]
          x_end_tleaf <- curve_data$tleaf[breadth_values_tleaf$idx_right]
          y_value_tleaf <- max(curve_data$pred_tleaf) * 0.95
        }
      } else if (model == "gam_tleaf_vpd") {
        # Fit GAM with tleaferature and VPD
        fit <- tryCatch(
          gam(photo ~ s(tleaf, bs='ts', k=5) + s(vpdl, k = 5, bs = 'ts'), 
              data = curve_data,
              select = TRUE,
              method = "REML"), 
          error = function(e) {
            print("BadFit: gam_tleaf_vpd")
            return(NULL)
          }
        )
        if (!is.null(fit)) {
          curve_data$pred_tleaf_vpdl <- predict(fit, newdata = curve_data)
          topt_values$topt_tleaf_vpdl <- gratia::smooth_estimates(fit, unconditional = TRUE, smooth = 's(tleaf)') %>% 
            filter(.estimate == max(.estimate)) %>% 
            pull(tleaf)
          breadth_values_tleaf_vpdl <- get_breadth(curve_data$pred_tleaf_vpdl, curve_data$tleaf, threshold=0.95)
          breadth_values$breadth_tleaf_vpdl <- breadth_values_tleaf_vpdl$breadth
          AIC_values$AIC_tleaf_vpdl <- AIC(fit)
          concurvity_values$concurvity_tleaf_vpdl <- concurvity(fit)
          
          # Calculate start and end points for the horizontal line
          x_start_tleaf_vpdl <- curve_data$tleaf[breadth_values_tleaf_vpdl$idx_left]
          x_end_tleaf_vpdl <- curve_data$tleaf[breadth_values_tleaf_vpdl$idx_right]
          y_value_tleaf_vpdl <- max(curve_data$pred_tleaf_vpdl) * 0.95
        }
      } else if (model == "gam_tleaf_vpd_cond") {
        # Fit GAM with tleaferature, VPD, and conductance
        fit <- tryCatch(
          gam(photo ~ s(tleaf, bs='ts', k=5) + s(vpdl, bs='ts', k=5) + s(cond, bs='ts', k=5), 
              data = curve_data,
              select = TRUE,
              method = "REML"), 
          error = function(e) {
            print("BadFit: gam_tleaf_vpd_cond")
            return(NULL)
          }
        )
        if (!is.null(fit)) {
          curve_data$pred_tleaf_vpdl_cond <- predict(fit, newdata = curve_data)
          topt_values$topt_tleaf_vpdl_cond <- gratia::smooth_estimates(fit, unconditional = TRUE, smooth = 's(tleaf)') %>% 
            filter(.estimate == max(.estimate)) %>% 
            pull(tleaf)
          breadth_values_tleaf_vpdl_cond <- get_breadth(curve_data$pred_tleaf_vpdl_cond, curve_data$tleaf, threshold=0.95)
          breadth_values$breadth_tleaf_vpdl_cond <- breadth_values_tleaf_vpdl_cond$breadth
          AIC_values$AIC_tleaf_vpdl_cond <- AIC(fit)
          concurvity_values$concurvity_tleaf_vpdl_cond <- concurvity(fit)
          
          # Calculate start and end points for the horizontal line
          x_start_tleaf_vpdl_cond <- curve_data$tleaf[breadth_values_tleaf_vpdl_cond$idx_left]
          x_end_tleaf_vpdl_cond <- curve_data$tleaf[breadth_values_tleaf_vpdl_cond$idx_right]
          y_value_tleaf_vpdl_cond <- max(curve_data$pred_tleaf_vpdl_cond) * 0.95
        }
      }
    }
    
    # Ensure the values are correctly extracted and assigned
    topt_tleaf <- ifelse(is.null(topt_values$topt_tleaf), NA, topt_values$topt_tleaf)
    topt_tleaf_cond <- ifelse(is.null(topt_values$topt_tleaf_cond), NA, topt_values$topt_tleaf_cond)
    topt_tleaf_vpdl <- ifelse(is.null(topt_values$topt_tleaf_vpdl), NA, topt_values$topt_tleaf_vpdl)
    topt_tleaf_vpdl_cond <- ifelse(is.null(topt_values$topt_tleaf_vpdl_cond), NA, topt_values$topt_tleaf_vpdl_cond)
    
    breadth_tleaf <- ifelse(is.null(breadth_values$breadth_tleaf), NA, breadth_values$breadth_tleaf)
    breadth_tleaf_cond <- ifelse(is.null(breadth_values$breadth_tleaf_cond), NA, breadth_values$breadth_tleaf_cond)
    breadth_tleaf_vpdl <- ifelse(is.null(breadth_values$breadth_tleaf_vpdl), NA, breadth_values$breadth_tleaf_vpdl)
    breadth_tleaf_vpdl_cond <- ifelse(is.null(breadth_values$breadth_tleaf_vpdl_cond), NA, breadth_values$breadth_tleaf_vpdl_cond)
    
    AIC_tleaf <- ifelse(is.null(AIC_values$AIC_tleaf), NA, AIC_values$AIC_tleaf)
    AIC_tleaf_cond <- ifelse(is.null(AIC_values$AIC_tleaf_cond), NA, AIC_values$AIC_tleaf_cond)
    AIC_tleaf_vpdl <- ifelse(is.null(AIC_values$AIC_tleaf_vpdl), NA, AIC_values$AIC_tleaf_vpdl)
    AIC_tleaf_vpdl_cond <- ifelse(is.null(AIC_values$AIC_tleaf_vpdl_cond), NA, AIC_values$AIC_tleaf_vpdl_cond)
    
    concurvity_tleaf <- ifelse(is.null(concurvity_values$concurvity_tleaf), NA, concurvity_values$concurvity_tleaf)
    concurvity_tleaf_cond <- ifelse(is.null(concurvity_values$concurvity_tleaf_cond), NA, concurvity_values$concurvity_tleaf_cond)
    concurvity_tleaf_vpdl <- ifelse(is.null(concurvity_values$concurvity_tleaf_vpdl), NA, concurvity_values$concurvity_tleaf_vpdl)
    concurvity_tleaf_vpdl_cond <- ifelse(is.null(concurvity_values$concurvity_tleaf_vpdl_cond), NA, concurvity_values$concurvity_tleaf_vpdl_cond)
    
    weib_topt <- ifelse(is.null(curve_data[[topt_col]]), NA, curve_data[[topt_col]])
    weib_breadth <- ifelse(is.null(curve_data[[breadth_col]]), NA, curve_data[[breadth_col]])
    
    results <- rbind(results, data.frame(curveid = curve_id, 
                                         topt_tleaf = topt_tleaf, 
                                         topt_tleaf_cond = topt_tleaf_cond, 
                                         topt_tleaf_vpdl = topt_tleaf_vpdl, 
                                         topt_tleaf_vpdl_cond = topt_tleaf_vpdl_cond, 
                                         topt_weib = weib_topt, 
                                         breadth_tleaf = breadth_tleaf, 
                                         breadth_tleaf_cond = breadth_tleaf_cond, 
                                         breadth_tleaf_vpdl = breadth_tleaf_vpdl, 
                                         breadth_tleaf_vpdl_cond = breadth_tleaf_vpdl_cond, 
                                         breadth_weib = weib_breadth, 
                                         AIC_tleaf = AIC_tleaf, 
                                         AIC_tleaf_cond = AIC_tleaf_cond, 
                                         AIC_tleaf_vpdl = AIC_tleaf_vpdl, 
                                         AIC_tleaf_vpdl_cond = AIC_tleaf_vpdl_cond, 
                                         concurvity_tleaf = max(concurvity_tleaf), 
                                         concurvity_tleaf_cond = max(concurvity_tleaf_cond), 
                                         concurvity_tleaf_vpdl = max(concurvity_tleaf_vpdl), 
                                         concurvity_tleaf_vpdl_cond = max(concurvity_tleaf_vpdl_cond)))
    
    # Plotting section
    p <- ggplot(curve_data, aes(x = tleaf, y = photo)) +
      geom_point() +
      geom_line(aes(y = pred_tleaf), color = "blue", linetype = "solid") +
      geom_line(aes(y = pred_tleaf_cond), color = "deeppink", linetype = "solid") +
      geom_line(aes(y = pred_tleaf_vpdl), color = "forestgreen", linetype = "solid") +
      geom_line(aes(y = pred_tleaf_vpdl_cond), color = "purple", linetype = "solid") +
      geom_vline(xintercept = weib_topt, linetype = "dashed", color = "grey") +
      geom_hline(yintercept = max(curve_data$photo, na.rm = TRUE) * 0.95, linetype = "dotted", color = "black") +
      theme_minimal() +
      labs(title = paste("Curve ID:", curve_id),
           x = "tleaferature",
           y = "Photosynthesis",
           color = "Model")
    
    # Add horizontal lines for breadth at 95% threshold
    if (!is.na(x_start_tleaf) & !is.na(x_end_tleaf)) {
      p <- p + geom_segment(aes(x = x_start_tleaf, y = y_value_tleaf, xend = x_end_tleaf, yend = y_value_tleaf), color = "blue")
    }
    if (!is.na(x_start_tleaf_cond) & !is.na(x_end_tleaf_cond)) {
      p <- p + geom_segment(aes(x = x_start_tleaf_cond, y = y_value_tleaf_cond, xend = x_end_tleaf_cond, yend = y_value_tleaf_cond), color = "red")
    }
    if (!is.na(x_start_tleaf_vpdl) & !is.na(x_end_tleaf_vpdl)) {
      p <- p + geom_segment(aes(x = x_start_tleaf_vpdl, y = y_value_tleaf_vpdl, xend = x_end_tleaf_vpdl, yend = y_value_tleaf_vpdl), color = "green")
    }
    if (!is.na(x_start_tleaf_vpdl_cond) & !is.na(x_end_tleaf_vpdl_cond)) {
      p <- p + geom_segment(aes(x = x_start_tleaf_vpdl_cond, y = y_value_tleaf_vpdl_cond, xend = x_end_tleaf_vpdl_cond, yend = y_value_tleaf_vpdl_cond), color = "purple")
    }
    
    print(p)
  }
  
  # Close the PDF file
  dev.off()
  
  # Return the results
  return(results)
}

# Example usage:
results <- fit_gams_and_plot(dat)
#filtered_results <- results %>%
#  rownames_to_column(var = "row_label")%>%
#  filter(str_detect(row_label, "worst"))

write.csv(results, "gam_fits_raw.csv")

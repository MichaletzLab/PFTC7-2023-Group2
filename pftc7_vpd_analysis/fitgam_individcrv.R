dat1 = dat%>%
  filter(curveid==c(3,1091))
# Load necessary libraries
library(mgcv)
library(ggplot2)
library(dplyr)
library(gratia)
library(data.table)

# Ensure your data.table is in place
dat1 <- as.data.table(dat1)

# Function to calculate breadth
get_breadth <- function(photo, tleaf, amax_fitted, threshold = 0.95) {
  # Find tleaf values around the peak
  idx_peak <- which.max(photo)
  idx_left <- max(1, idx_peak - 1)
  idx_right <- min(length(photo), idx_peak + 1)
  
  # Calculate breadth
  breadth_left <- tleaf[idx_peak] - tleaf[idx_left]
  breadth_right <- tleaf[idx_right] - tleaf[idx_peak]
  breadth <- breadth_left + breadth_right
  
  # Scale breadth by the fitted Amax value
  breadth_scaled <- breadth * (threshold * amax_fitted)
  
  return(breadth_scaled)
}

# Function to fit GAMs and plot each curve individually
fit_gams_and_plot <- function(data, curveid_col = "curveid", x_col = "tleaf", y_col = "photo", cond_col = "cond", topt_col = "weib.topt", breadth_col = "weib.breadth") {
  models <- c("gam_temp_cond", "gam_temp")
  results <- data.frame(curveid = numeric(), topt_temp = numeric(), topt_temp_cond = numeric(), topt_weib = numeric(), breadth_temp = numeric(), breadth_temp_cond = numeric(), breadth_weib = numeric())
  
  # Create a PDF file to save the plots
  pdf("gam_fits.pdf")
  
  for (curve_id in unique(data[[curveid_col]])) {
    curve_data <- subset(data, data[[curveid_col]] == curve_id) # Subset data for the current curveID
    print(paste("Curve ID:", curve_id, ", Number of observations:", nrow(curve_data)))  # Debug
    
    # Initialize prediction columns
    curve_data$pred_temp <- NA
    curve_data$pred_temp_cond <- NA
    
    topt_values <- list(topt_temp = NA, topt_temp_cond = NA)
    breadth_values <- list(breadth_temp = NA, breadth_temp_cond = NA)
    
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
          amax_fitted <- fitted(fit, type = "lpmatrix")[1]  # Fitted Amax value
          breadth_values$breadth_temp_cond <- get_breadth(curve_data$photo, curve_data$tleaf, amax_fitted)
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
          amax_fitted <- fitted(fit, type = "lpmatrix")[1]  # Fitted Amax value
          breadth_values$breadth_temp <- get_breadth(curve_data$photo, curve_data$tleaf, amax_fitted)
        }
      }
    }
    
    topt_weib <- mean(curve_data[[topt_col]], na.rm = TRUE)
    topt_values$topt_weib <- topt_weib
    breadth_weib <- mean(curve_data[[breadth_col]])
    breadth_values$breadth_weib <- breadth_weib
    
    # Save the Topt and breadth values in the results dataframe
    results <- rbind(results, data.frame(curveid = curve_id, 
                                         topt_temp = topt_values$topt_temp, 
                                         topt_temp_cond = topt_values$topt_temp_cond, 
                                         topt_weib = topt_weib,
                                         breadth_temp = breadth_values$breadth_temp,
                                         breadth_temp_cond = breadth_values$breadth_temp_cond,
                                         breadth_weib = breadth_weib))
    
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
    
    # Print plot to PDF
    print(plot)
  }
  
  # Close the PDF device
  dev.off()
  
  return(results)
}

# Call the function with your data
results <- fit_gams_and_plot(dat1)

# Print the results dataframe
print(results)

# Confirm the PDF creation
if (file.exists("gam_fits.pdf")) {
  message("PDF successfully created.")
} else {
  message("PDF creation failed. Please check the code and directory path.")
}

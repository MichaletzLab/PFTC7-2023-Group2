cut.hooks <- function(data, pdf_file) {
  data$no_modes <- NA
  
  final_data <- data.frame()  # Initialize empty dataframe to store final data
  
  # Open a PDF device
  pdf(pdf_file)
  
  for(i in unique(data$curveID)) {
    d1 <- subset(data, curveID == i)  # Subset data for the current curve
    
    # Filter data for Tleaf > 25
    d1_high <- subset(d1, Tleaf > 28)
    
    # Check if there are enough data points above the threshold
    if (nrow(d1_high) < 10) {  # Adjust threshold as needed
      next  # Skip this curve if there are not enough data points
    }
    
    # Fit segmented regression model with multiple breakpoints
    segmented_model <- segmented(lm(A ~ Tleaf, data = d1_high), seg.Z = ~ Tleaf, nk = 2)  # You can adjust nk as needed
    
    # Extract breakpoints
    breakpoints <- segmented_model$psi[, "Est."]
    
    # Print the breakpoints for the current curve
    cat("Breakpoints for Curve ID", i, ":", breakpoints, "\n")
    
    # Check if a breakpoint is found beyond Tleaf 25
    if (length(breakpoints) > 0 && any(breakpoints > 28)) {
      # Get the last breakpoint beyond Tleaf 25
      last_breakpoint <- max(breakpoints[breakpoints > 28])
      
      # Find the corresponding Tleaf value for the last breakpoint
      last_breakpoint_Tleaf <- d1$Tleaf[which.max(d1$Tleaf > last_breakpoint)]
      
      # Truncate the data after the point of the last breakpoint beyond Tleaf 25
      truncated_data <- d1_high[d1_high$Tleaf <= last_breakpoint_Tleaf, ]
      
      # Plot the truncated curve with segmented lines
      plot(d1$Tleaf, d1$A, xlab = "Tleaf", ylab = "A", main = paste("Truncated Curve ID:", i))
      lines(d1_high$Tleaf, d1_high$A, col = "red", lwd = 2)
      abline(v = last_breakpoint_Tleaf, col = "blue", lty = 2)
      
      # Add truncated data to final_data
      final_data <- bind_rows(final_data, truncated_data)
    } else {
      # If no breakpoint is found beyond Tleaf 25, keep the original data
      final_data <- bind_rows(final_data, d1_high)
      
      # Plot the original curve with segmented lines
      plot(d1$Tleaf, d1$A, xlab = "Tleaf", ylab = "A", main = paste("Original Curve ID:", i))
      lines(d1_high$Tleaf, predict(segmented_model), col = "red", lwd = 2)
    }
    
    # Combine the truncated data with the original data for Tleaf <= 25
    final_data <- bind_rows(final_data, subset(d1, Tleaf <= 28))
    
    # Check if a breakpoint is found
    if (length(breakpoints) > 0) {
      # Discard curves where a breakpoint is found
      data[data$curveID == i,]$no_modes <- 0
    } else {
      data[data$curveID == i,]$no_modes <- 1
    }
  }
  
  # Close the PDF device
  dev.off()
  
  return(final_data)
}


#at.subset3 = cut.hooks(at.subset2,"truncate_plots.pdf")
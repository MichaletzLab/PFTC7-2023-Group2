discard.hooks = function(data, threshold = 0) { 
  data$no_modes = NA
  curves_to_discard <- c() 
  
  for(curve_id in unique(data$curveID)) {
    
    d1 = subset(data, curveID == curve_id)
    
    # Subset the last 20 data points
    last_20_points = tail(d1, 50)
    
    # Calculate the slope between the first and last of the last 20 points
    slope_last_20 = (last_20_points$A[nrow(last_20_points)] - last_20_points$A[1]) / (last_20_points$Tleaf[nrow(last_20_points)] - last_20_points$Tleaf[1])
    
    # Check if the slope curvature exceeds the threshold
    if (slope_last_20 > threshold) {
      curves_to_discard <- c(curves_to_discard, curve_id)
    }
    
    # Plot the curve with last 20 points highlighted
    plot(d1$Tleaf, d1$A, type = "l", main = paste("Curve ID:", curve_id), xlab = "Tleaf", ylab = "A")
    points(last_20_points$Tleaf, last_20_points$A, col = "blue", pch = 16)
    
    # Mark the end of the curve if it's being discarded
    if (curve_id %in% curves_to_discard) {
      points(last_20_points$Tleaf[nrow(last_20_points)], last_20_points$A[nrow(last_20_points)], col = "red", pch = 16)
    }
    
    # Add a legend
    legend("bottomleft", legend = c("Data", "Last 20 Points", "End of Curve"), col = c("black", "blue", "red"), pch = c(NA, 16, 16), lty = c(1, NA, NA))
    
    # Add a pause to see each plot (optional)
    Sys.sleep(1)
    
  }
  
  print(paste("Number of curves discarded due to curvature at the end:", length(curves_to_discard)))
  
  # Exclude curves with significant curvature at the end
  data_filtered = data[!data$curveID %in% curves_to_discard, ]
  
  return(data_filtered)
}




#Test
#test = discard.hooks(test123,0.01) #0.01 works well enough
#unique(data_cut$curveID)





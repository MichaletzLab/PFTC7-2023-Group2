# Function to remove multimodal curves from dataset
cut.multimodal = function(data) {
  
  data$no_modes <- NA
  
  for (i in unique(data$curveID)) {
    d1 <- subset(data, curveID == i)
    
    # Check for NA, Inf, or -Inf values and remove them
    #d1 <- d1[complete.cases(d1), ]
    #d1 <- d1[is.finite(d1$Tleaf) & is.finite(d1$A), ]
    
    # Check if there are enough data points
    if (nrow(d1) < 3) {
      data[data$curveID == i, ]$no_modes <- NA
      next
    }
    
    # Perform spline smoothing
    tryCatch({
      z <- smooth.spline(x = d1$Tleaf, y = d1$A, nknots = 10)
      
      Sx <- seq(min(d1$Tleaf) + 1, max(d1$Tleaf) - 1, 0.01)
      len <- length(Sx)
      pred.spline <- predict(z, Sx)
      d0 <- data.frame(pred.spline)
      pred.prime <- predict(z, Sx, deriv = 1)
      d0$local_max <- c(pred.prime$y[1:(len - 1)] > 0 & pred.prime$y[2:len] < 0, FALSE)
      
      Mins <- which(d0$local_max == TRUE)
      data[data$curveID == i, ]$no_modes <- length(Mins)
      print(paste("CurveID:", i, "- Spline fitting successful."))
      
    }, error = function(e) {
      # Handle error by setting no_modes to NA
      data[data$curveID == i, ]$no_modes <- NA
    })
  }
  
  data.cut <- subset(data, no_modes == 1)
  return(data.cut)
}
# Function to remove data w/ poor r-squared or with fitted T_opt
# too close to Tmin or Tmax
cut.topt.out.of.bounds = function(data) {
  
  data2 = data %>% dplyr::select(Tleaf, A, curveID)
  
  results = fit_weibull_breadth(data2 %>% dplyr::select(Tleaf, A, curveID), x = "Tleaf", y = "A")
  
  # For each curve, check Topt > Tmin and Topt < Tmax
  results$keep = NA
  for (i in unique(data$curveID)) {
    curdat = subset(data, curveID == i)
    Tmin = min(curdat$Tleaf)
    Tmax = max(curdat$Tleaf)
    j = which(results$curveID == i)
    Topt = results[j,]$T_opt
    r_sq = results[j,]$r_sq
    
    results[j,]$keep = Topt > Tmin+3 & Topt < Tmax-3 & r_sq > 0.9
  }
  
  
  results = subset(results, keep == T)
  
  data = subset(data, curveID %in% results$curveID)
  return(data)
  
}


#grrr=cut.topt.out.of.bounds(at.subset)
#
# at.subset = at.corr
#
#at.subset1 = cut.multimodal(at.subset)
#unique(at.subset1$curveID)
#at.subset2 = cut.topt.out.of.bounds(at.subset1)
#unique(at.subset2$curveID)


#results <- fit_weibull(at.subset %>% select(Tleaf, A, curveID), x = "Tleaf", y = "A")#
# meta = unique(bind_cols(curveID = as.numeric(at.subset$curveID),at.subset[,284:309]))
#
# results_meta = merge(results, meta, by="curveID")

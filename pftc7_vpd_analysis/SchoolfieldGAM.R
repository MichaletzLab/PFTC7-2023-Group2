dat1 <- dat%>%
  filter(curveid==c(3,6,9,1091))

schoolfield.fit = read.csv("outputs/discard.schoolfield.SANW.csv")%>%
  rename(curveid = curveID)%>%
  select(curveid, T_opt,J_ref, E, E_D,breadth)%>%
  rename(school_breadth = breadth)%>%
  rename(T_opt.s = T_opt)
mjcschoolfield.fit = read.csv("MJCschoolfield.fit.csv")%>%
  rename(curveid = curveID)%>%
  select(curveid, T_opt,J_ref, E, E_D,breadth)%>%
  rename(mjcschool_breadth = breadth)%>%
  rename(T_opt.m = T_opt)
weibull.fit = read.csv("outputs/discard.weibull.SANW.csv")%>%
  rename(curveid = curveID)%>%
  rename(T_opt.w = T_opt)

# Define the Schoolfield model function
schoolfield <- function(temp, J_ref, E, E_D, T_opt, T_ref = 25) {
  k <- 8.62e-5            # Boltzmann's constant (eV K-1)
  temp <- temp + 273.15    # Convert to Kelvin
  T_ref <- T_ref + 273.15  # Convert to Kelvin
  T_opt <- T_opt + 273.15  # Convert to Kelvin
  return(J_ref * exp(E * (1/(k*T_ref) - 1/(k*temp))) / (1 + E/(E_D - E) * exp((E_D/k) * (1/T_opt - 1/temp))))
}
mjcschoolfield <- function(temp, J_ref, E, E_D, T_opt, T_ref = 25, Tair) {
  k <- 8.62e-5            # Boltzmann's constant (eV K-1)
  temp <- temp + 273.15    # Convert to Kelvin
  T_ref <- T_ref + 273.15  # Convert to Kelvin
  T_opt <- T_opt + 273.15  # Convert to Kelvin
  Tair <- Tair + 273.15    # Convert Tair to Kelvin
  a <- 610.78
  b <- 237.3
  c <- 17.2694
  d <- 0.35
  VPD <- (a * exp(temp / ((temp + b) * c))) - (d * a * exp(Tair / ((Tair + b) * c)))
  return((J_ref * exp(E * (1 / (k * T_ref) - 1 / (k * temp))) / (1 + E / (E_D - E) * exp((E_D / k) * (1 / T_opt - 1 / temp)))) * VPD)
}

# Merge dat1 with schoolfield.fits based on curveid
dats <- dat %>%
  left_join(schoolfield.fit, by = "curveid")%>%
  left_join(weibull.fit, by = "curveid")

# Calculate the predicted values using the Schoolfield model
dats <- dats %>%
  mutate(predicted = schoolfield(tleaf, J_ref, E, E_D, T_opt.s))
# Calculate the predicted values using the MJCschoolfield model
dats <- dats %>%
  mutate(predicted.mjc = mjcschoolfield(tleaf, J_ref, E, E_D, T_opt.w,Tair = tair))

# Define the function to fit the GAM model and extract T_opt and curve breadth
fit_gam_with_schoolfield <- function(data, x_var, y_var, curveID_var) {
  results <- data.frame()
  pdf("outputs/gam_fit_compare.pdf")
  
  for (curve_id in unique(data[[curveID_var]])) {
    print(paste("Fitting curve ID:", curve_id))
    curve_data <- data[data[[curveID_var]] == curve_id, ]  # Subset using data frame indexing
    
    if (nrow(curve_data) < 4) {
      message(paste("Skipping curve ID:", curve_id, "because it has fewer than 4 observations."))
      next  # Skip to the next curve if there are fewer than 4 observations
    }
    
    # Fit the GAM models
    formula0 <- as.formula(paste(y_var, "~ s(tleaf) + predicted"))
    formula1 <- as.formula(paste(y_var, "~ s(cond) + predicted"))
    formula2 <- as.formula(paste(y_var, "~ s(cond) + s(vpdl) + predicted"))
    formula3 <- as.formula(paste(y_var, "~ s(tleaf) + predicted.mjc"))
    formula4 <- as.formula(paste(y_var,"~s(tleaf)"))
    formula5 <- as.formula(paste(y_var,"~s(tleaf)+s(vpdl)"))
    formula6 <- as.formula(paste(y_var,"~s(tleaf)+s(cond)"))
    
    tryCatch({
      #####################
      #####################
      # Fit the GAM model with s(tleaf)
      gam_fit0 <- gam(formula0, data = curve_data, method = "REML")
      conc0 <- concurvity(gam_fit0, full = TRUE)
      
      # Extract fitted values from GAM model with s(tleaf)
      fitted_values0 <- predict(gam_fit0, newdata = curve_data, type = "response")
      max_photo_gam0 <- max(fitted_values0)
      
      # Extract new T_opt for GAM model with s(tleaf)
      new_T_opt0 <- curve_data[[x_var]][which.max(fitted_values0)]
      
      # Calculate 85% of max photo for GAM model with s(tleaf)
      threshold_85_gam0 <- 0.85 * max_photo_gam0
      within_85_gam0 <- curve_data[[x_var]][fitted_values0 >= threshold_85_gam0]
      curve_breadth_gam0 <- max(within_85_gam0) - min(within_85_gam0)
      #####################
      #####################
      # Fit the GAM model with s(cond)
      gam_fit1 <- gam(formula1, data = curve_data, method = "REML")
      conc1 <- concurvity(gam_fit1, full = TRUE)
      
      # Extract fitted values from GAM model with s(cond)
      fitted_values1 <- predict(gam_fit1, newdata = curve_data, type = "response")
      max_photo_gam1 <- max(fitted_values1)
      
      # Extract new T_opt for GAM model with s(cond)
      new_T_opt1 <- curve_data[[x_var]][which.max(fitted_values1)]
      
      # Calculate 85% of max photo for GAM model with s(cond)
      threshold_85_gam1 <- 0.85 * max_photo_gam1
      within_85_gam1 <- curve_data[[x_var]][fitted_values1 >= threshold_85_gam1]
      curve_breadth_gam1 <- max(within_85_gam1) - min(within_85_gam1)
      #####################
      #####################
      # Fit the GAM model with s(cond) + s(vpdl)
      gam_fit2 <- gam(formula2, data = curve_data, method = "REML")
      conc2 <- concurvity(gam_fit2, full = TRUE)
      
      # Extract fitted values from GAM model with s(cond) + s(vpdl)
      fitted_values2 <- predict(gam_fit2, newdata = curve_data, type = "response")
      max_photo_gam2 <- max(fitted_values2)
      
      # Extract new T_opt for GAM model with s(cond) + s(vpdl)
      new_T_opt2 <- curve_data[[x_var]][which.max(fitted_values2)]
      
      # Calculate 85% of max photo for GAM model with s(cond) + s(vpdl)
      threshold_85_gam2 <- 0.85 * max_photo_gam2
      within_85_gam2 <- curve_data[[x_var]][fitted_values2 >= threshold_85_gam2]
      curve_breadth_gam2 <- max(within_85_gam2) - min(within_85_gam2)
      #####################
      #####################
      gam_fit3 <- gam(formula3, data = curve_data, method = "REML")
      conc3 <- concurvity(gam_fit3, full = TRUE)
      
      # Extract fitted values from GAM model with s(tleaf)
      fitted_values3 <- predict(gam_fit3, newdata = curve_data, type = "response")
      max_photo_gam3 <- max(fitted_values3)
      
      # Extract new T_opt for GAM model with s(tleaf)
      new_T_opt3 <- curve_data[[x_var]][which.max(fitted_values3)]
      
      # Calculate 85% of max photo for GAM model with s(tleaf)
      threshold_85_gam3 <- 0.85 * max_photo_gam3
      within_85_gam3 <- curve_data[[x_var]][fitted_values3 >= threshold_85_gam3]
      curve_breadth_gam3 <- max(within_85_gam3) - min(within_85_gam3)
      #####################
      #####################
      gam_fit4 <- gam(formula4, data = curve_data, method = "REML")
      conc4 <- concurvity(gam_fit4, full = TRUE)
      
      # Extract fitted values from GAM model with s(tleaf)
      fitted_values4 <- predict(gam_fit4, newdata = curve_data, type = "response")
      max_photo_gam4 <- max(fitted_values4)
      
      # Extract new T_opt for GAM model with s(tleaf)
      new_T_opt4 <- curve_data[[x_var]][which.max(fitted_values4)]
      
      # Calculate 85% of max photo for GAM model with s(tleaf)
      threshold_85_gam4 <- 0.85 * max_photo_gam4
      within_85_gam4 <- curve_data[[x_var]][fitted_values4 >= threshold_85_gam4]
      curve_breadth_gam4 <- max(within_85_gam4) - min(within_85_gam4)
      #####################
      #####################
      gam_fit5 <- gam(formula5, data = curve_data, method = "REML")
      conc5 <- concurvity(gam_fit5, full = TRUE)
      
      # Extract fitted values from GAM model with s(tleaf)
      fitted_values5 <- predict(gam_fit5, newdata = curve_data, type = "response")
      max_photo_gam5 <- max(fitted_values5)
      
      # Extract new T_opt for GAM model with s(tleaf)
      new_T_opt5 <- curve_data[[x_var]][which.max(fitted_values5)]
      
      # Calculate 85% of max photo for GAM model with s(tleaf)
      threshold_85_gam5 <- 0.85 * max_photo_gam5
      within_85_gam5 <- curve_data[[x_var]][fitted_values5 >= threshold_85_gam5]
      curve_breadth_gam5 <- max(within_85_gam5) - min(within_85_gam5)
      #####################
      #####################
      gam_fit6 <- gam(formula6, data = curve_data, method = "REML")
      conc6 <- concurvity(gam_fit6, full = TRUE)
      
      # Extract fitted values from GAM model with s(tleaf)
      fitted_values6 <- predict(gam_fit6, newdata = curve_data, type = "response")
      max_photo_gam6 <- max(fitted_values6)
      
      # Extract new T_opt for GAM model with s(tleaf)
      new_T_opt6 <- curve_data[[x_var]][which.max(fitted_values6)]
      
      # Calculate 85% of max photo for GAM model with s(tleaf)
      threshold_85_gam6 <- 0.85 * max_photo_gam6
      within_85_gam6 <- curve_data[[x_var]][fitted_values6 >= threshold_85_gam6]
      curve_breadth_gam6 <- max(within_85_gam6) - min(within_85_gam6)
      #####################
      #####################
      
      # Calculate 85% of max photo for Schoolfield model
      predicted_values <- curve_data$predicted
      max_photo_schoolfield <- max(predicted_values)
      threshold_85_schoolfield <- 0.85 * max_photo_schoolfield
      within_85_schoolfield <- curve_data[[x_var]][predicted_values >= threshold_85_schoolfield]
      curve_breadth_schoolfield <- max(within_85_schoolfield) - min(within_85_schoolfield)
      
      predicted.mjc_values <- curve_data$predicted.mjc
      max_photo_schoolfield.mjc <- max(predicted.mjc_values)
      threshold_85_schoolfield.mjc <- 0.85 * max_photo_schoolfield.mjc
      within_85_schoolfield.mjc <- curve_data[[x_var]][predicted.mjc_values >= threshold_85_schoolfield.mjc]
      curve_breadth_schoolfield.mjc <- max(within_85_schoolfield.mjc) - min(within_85_schoolfield.mjc)
      
      # Extract T_opt_schoolfield from left-joined data (assuming it's named T_opt_schoolfield)
      T_opt_schoolfield <- mean(curve_data$T_opt.s)
      T_opt_schoolfield.mjc <- mean(curve_data$T_opt.m)
      T_opt_weibull <- mean(curve_data$T_opt.w)
      # Plotting
      plot <- ggplot(curve_data, aes_string(x = x_var, y = y_var)) +
        geom_point() +
        geom_smooth(aes(color = "purple4"), method = "gam", se = TRUE, formula = formula4) +
        geom_smooth(aes(color = "turquoise4"), method = "gam", se = TRUE, formula = formula5) +
        geom_smooth(aes(color = "lightblue3"), method = "gam", se = TRUE, formula = formula6) +
        geom_smooth(aes(color = "darkgreen"), method = "gam", se = TRUE, formula = formula0) +
        geom_smooth(aes(color = "blue"), method = "gam", se = TRUE, formula = formula1) +
        geom_smooth(aes(color = "purple"), method = "gam", se = TRUE, formula = formula2) +
        geom_smooth(aes(color = "turquoise"), method = "gam", se = TRUE, formula = formula3) +
        geom_line(aes(color = "purple4", y = fitted_values4)) +
        geom_line(aes(color = "turquoise4", y = fitted_values5)) +
        geom_line(aes(color = "lightblue3", y = fitted_values6)) +
        geom_line(aes(color = "darkgreen", y = fitted_values0)) +
        geom_line(aes(color = "blue", y = fitted_values1)) +
        geom_line(aes(color = "purple", y = fitted_values2)) +
        geom_line(aes(color = "turquoise", y = fitted_values3)) +
        geom_line(aes(color = "red", y = predicted_values)) +
        #geom_line(aes(color = "maroon", y = predicted.mjc_values)) +
        geom_vline(xintercept = new_T_opt4, linetype = "dashed", color = "purple4") +
        geom_vline(xintercept = new_T_opt5, linetype = "dashed", color = "turquoise4") +
        geom_vline(xintercept = new_T_opt6, linetype = "dashed", color = "lightblue3") +
        geom_vline(xintercept = new_T_opt0, linetype = "dashed", color = "darkgreen") +
        geom_vline(xintercept = new_T_opt1, linetype = "dashed", color = "blue") +
        geom_vline(xintercept = new_T_opt2, linetype = "dashed", color = "purple") +
        geom_vline(xintercept = new_T_opt3, linetype = "dashed", color = "turquoise") +
        geom_point(aes(x = new_T_opt4, y = threshold_85_gam4, color = "purple4"), size = 3, shape = 16) +
        geom_point(aes(x = new_T_opt5, y = threshold_85_gam5, color = "turquoise4"), size = 3, shape = 16) +
        geom_point(aes(x = new_T_opt6, y = threshold_85_gam6, color = "lightblue3"), size = 3, shape = 16) +
        geom_point(aes(x = new_T_opt0, y = threshold_85_gam0, color = "darkgreen"), size = 3, shape = 16) +
        geom_point(aes(x = new_T_opt1, y = threshold_85_gam1, color = "blue"), size = 3, shape = 16) +
        geom_point(aes(x = new_T_opt2, y = threshold_85_gam2, color = "purple"), size = 3, shape = 16) +
        geom_point(aes(x = new_T_opt3, y = threshold_85_gam3, color = "turquoise"), size = 3, shape = 16) +
        geom_point(aes(x = T_opt_schoolfield, y = threshold_85_schoolfield, color = "red"), size = 3, shape = 16) +
        geom_vline(xintercept = T_opt_schoolfield, linetype = "dashed", color = "red") +
        geom_vline(xintercept = T_opt_schoolfield.mjc, linetype = "dashed", color = "maroon") +
        geom_vline(xintercept = T_opt_weibull, linetype = "dashed", color = "coral") +
        ggtitle(paste("Curve ID:", curve_id)) +
        theme_minimal() +
        labs(color = "Model") +  # Legend title
        scale_color_manual(
          values = c(
            "darkgreen" = "darkgreen", 
            "blue" = "blue", 
            "purple" = "purple", 
            "turquoise" = "turquoise",
            "purple4" = "purple4", 
            "turquoise4" = "turquoise4", 
            "lightblue3" = "lightblue3", 
            "red" = "red",
            "maroon" = "maroon", 
            "coral" = "coral"
          ),
          breaks = c(
            "purple4", 
            "turquoise4", 
            "lightblue3",
            "darkgreen", 
            "blue", 
            "purple", 
            "turquoise", 
            "red", 
            "maroon", 
            "coral"
          ),
          labels = c(
            "GAM with s(tleaf)", 
            "GAM with s(tleaf) + s(gsw)", 
            "GAM with s(tleaf) + s(vpd)", 
            "SS GAM with s(tleaf)", 
            "SS GAM with s(gsw)", 
            "SS GAM with s(gsw) + s(vpd)", 
            "SS x vpd GAM with s(tleaf)", 
            "Schoolfield", 
            "SS x vpd", 
            "Weibull"
          )
        )  # Legend labels
      print(plot)
      
      # Collect results
      result_row <- data.frame(
        curve_id = curve_id,
        schoolfield.T_opt = T_opt_schoolfield,
        schoolfield.mjc.T_opt = T_opt_schoolfield.mjc,
        weibull.T_opt = T_opt_weibull,
        new_T_opt0 = new_T_opt0,
        new_T_opt1 = new_T_opt1,
        new_T_opt2 = new_T_opt2,
        new_T_opt3 = new_T_opt3,
        new_T_opt4 = new_T_opt4,
        new_T_opt5 = new_T_opt5,
        new_T_opt6 = new_T_opt6,
        curve_breadth_gam0 = curve_breadth_gam0,
        curve_breadth_gam1 = curve_breadth_gam1,
        curve_breadth_gam2 = curve_breadth_gam2,
        curve_breadth_gam3 = curve_breadth_gam3,
        curve_breadth_gam4 = curve_breadth_gam4,
        curve_breadth_gam5 = curve_breadth_gam5,
        curve_breadth_gam6 = curve_breadth_gam6,
        curve_breadth_schoolfield = curve_breadth_schoolfield,
        curve_breadth_schoolfield.mjc = curve_breadth_schoolfield.mjc,
        AIC0 = AIC(gam_fit0),
        AIC1 = AIC(gam_fit1),
        AIC2 = AIC(gam_fit2),
        AIC3 = AIC(gam_fit3),
        AIC4 = AIC(gam_fit4),
        AIC5 = AIC(gam_fit5),
        AIC6 = AIC(gam_fit6),
        rsq0 = summary(gam_fit0)$r.sq,
        rsq1 = summary(gam_fit1)$r.sq,
        rsq2 = summary(gam_fit2)$r.sq,
        rsq3 = summary(gam_fit3)$r.sq,
        rsq4 = summary(gam_fit4)$r.sq,
        rsq5 = summary(gam_fit5)$r.sq,
        rsq6 = summary(gam_fit6)$r.sq,
        concurvity_full0 = conc0,
        concurvity_full1 = conc1,
        concurvity_full2 = conc2,
        concurvity_full3 = conc3,
        concurvity_full4 = conc4,
        concurvity_full5 = conc5,
        concurvity_full6 = conc6
      )
      
      results <- rbind(results, result_row)
      
    }, error = function(e) {
      message(paste("Failed to fit curve ID:", curve_id, "with error:", e$message))
    })
  }
  dev.off()
  return(results)
}


results <- fit_gam_with_schoolfield(dats, "tleaf", "photo", "curveid")
results <- results%>%
  rename(curveid = curve_id)%>%
  left_join(schoolfield.fit, by = "curveid")

results <- results %>% as_tibble()

results_filtered <- results %>%
  group_by(curveid) %>%
  filter(concurvity_full1.s.cond. == max(concurvity_full1.s.cond.))

write.csv(results_filtered, "GAM.compare.csv")
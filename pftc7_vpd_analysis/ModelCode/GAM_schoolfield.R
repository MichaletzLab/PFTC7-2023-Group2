################################################################
# Schoolfield Only GAM:::: 
################################################################

# Custom smooth constructor for the Schoolfield model
schoolfield_smooth <- function(temp, J_ref, E, E_D, T_opt, k = 8.62e-5) {
  temp <- temp + 273.15  # Convert to Kelvin
  T_ref <- 25 + 273.15  # Reference temperature (Kelvin)
  T_opt <- T_opt + 273.15  # Optimal temperature (Kelvin)
  return(J_ref * exp(E * (1/(k*T_ref) - 1/(k*temp))) / 
           (1 + E/(E_D - E) * exp((E_D/k) * (1/T_opt - 1/temp))))
}

# Define a function to generate the Schoolfield basis
schoolfield_basis <- function(data, x_var, J_ref, E, E_D, T_opt) {
  x <- data[[x_var]]
  basis <- schoolfield_smooth(x, J_ref, E, E_D, T_opt)  # Ensure correct function call here
  return(basis)
}

# Wrapper function to fit the GAM model with Schoolfield parameters
fit_schoolfieldONLY_gam <- function(data, y_var, x_var, T_ref = 25, start_params) {
  # Objective function to minimize (negative log-likelihood)
  obj_fun <- function(params) {
    J_ref <- params[1]
    E <- params[2]
    E_D <- params[3]
    T_opt <- params[4]
    data$schoolfield_base <- schoolfield_basis(data, x_var, J_ref, E, E_D, T_opt, T_ref)
    gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
    return(-logLik(gam_fit))
  }
  
  # Optimize the parameters
  opt_res <- optim(start_params, obj_fun, method = "L-BFGS-B", lower = c(0, 0, 0.2, 0), upper = c(40, 5, 15, 40))
  
  # Extract the optimized parameters
  opt_params <- opt_res$par
  J_ref <- opt_params[1]
  E <- opt_params[2]
  E_D <- opt_params[3]
  T_opt <- opt_params[4]
  
  # Print the optimized parameters
  cat("Optimized J_ref:", J_ref, "\n")
  cat("Optimized E:", E, "\n")
  cat("Optimized E_D:", E_D, "\n")
  cat("Optimized T_opt:", T_opt, "\n")
  
  # Fit the GAM model with optimized parameters
  data$schoolfield_base <- schoolfield_basis(data, x_var, J_ref, E, E_D, T_opt, T_ref)
  gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
  
  return(list(gam_fit = gam_fit, opt_params = opt_params))
}


#start_params <- c(J_ref = 5, E = 0.5, E_D = 1.5, T_opt = 20)  # These were determined by the average of the nls.multstart results
#data <- dat1 %>% select(tleaf, photo)
#fit_results_only_schoolfield <- fit_schoolfieldONLY_gam(data, "photo", "tleaf", T_ref = 25, start_params = start_params)

#summary(fit_results_only_schoolfield$gam_fit)




####################Partial Effects Plot:::::::::::::::
# Define a range of tleaf values for prediction
# tleaf_range <- seq(min(data$tleaf), max(data$tleaf), length.out = 100)
# 
# # Create a new dataframe for prediction
# prediction_data <- data.frame(tleaf = tleaf_range)
# 
# # Add the schoolfield_basis to the prediction data using Celsius values
# prediction_data$schoolfield_base <- schoolfield_smooth(prediction_data$tleaf, 
#                                                        fit_results_only_schoolfield$opt_params[1], 
#                                                        fit_results_only_schoolfield$opt_params[2], 
#                                                        fit_results_only_schoolfield$opt_params[3], 
#                                                        fit_results_only_schoolfield$opt_params[4],
#                                                        T_ref = 25)  # Ensure T_ref is in Celsius

# # Predict the responses and standard errors
# predictions <- predict(fit_results_only_schoolfield$gam_fit, 
#                        newdata = prediction_data, type = "response", se.fit = TRUE)

# Add the predictions and confidence intervals to the prediction data
# prediction_data$predicted <- predictions$fit
# prediction_data$predicted_upper <- predictions$fit + 1.96 * predictions$se.fit
# prediction_data$predicted_lower <- predictions$fit - 1.96 * predictions$se.fit

# Plot the partial effects with error bands
# p <- ggplot(prediction_data, aes(x = tleaf, y = predicted)) +
#   geom_ribbon(aes(ymin = predicted_lower, ymax = predicted_upper), alpha = 0.2, fill = "blue") +
#   geom_line(color = "blue") +
#   labs(x = "Leaf temperature (°C)", 
#        y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m^-2, s^-1, ")")),
#        title = "Partial Effects Plot with Confidence Intervals") +
#   theme_minimal()


# Save the plot
# ggsave(p, 
#        filename = paste0("pftc7_vpd_analysis/figures/schoolfieldONLY_GAM_", Sys.Date(), ".png"),
#        device = "png",
#        width = 8,
#        height = 6,
#        units = 'in',
#        dpi = 300)


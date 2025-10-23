#####################
#Schoolfield in gam with cond and tleaf
#####################

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

# Wrapper function to fit the GAM model with Schoolfield parameters and smoothing term for conductance
fit_schoolfield.tleaf.cond_gam <- function(data, x_var, y_var, cond_var,  T_ref = 25, start_params) {
  # Objective function to minimize (negative log-likelihood)
  obj_fun <- function(params) {
    J_ref <- params[1]
    E <- params[2]
    E_D <- params[3]
    T_opt <- params[4]
    data$schoolfield_basis <- schoolfield_basis(data, x_var, J_ref, E, E_D, T_opt)
    gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_basis + s(", cond_var, ", k=5) + s(", x_var, ", k=5)")), data = data, method = "REML")
    return(-logLik(gam_fit))
  }
  
  # Optimize the parameters
  opt_res <- optim(start_params, obj_fun, method = "L-BFGS-B", lower = c(0, 0, 0.2, 0), upper = c(20, 2, 15, 50))
  
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
  data$schoolfield_basis <- schoolfield_basis(data, x_var, J_ref, E, E_D, T_opt)
  gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_basis + s(", cond_var, ", k=5) +s(", x_var, ", k=5)")), data = data, method = "REML")
  model_aic <- AIC(gam_fit)
  
  # Check concurvity
  model_concurvity <- concurvity(gam_fit)
  
  return(list(gam_fit = gam_fit, opt_params = opt_params, AIC = model_aic, concurvity = model_concurvity))
}


#start_params <- c(J_ref = 35, E = 0.5, E_D = 1.5, T_opt = 20) #these were determined by the average of the nls.multstart results
#data <- dat %>% select(tleaf, photo, cond)
#GAM.SS.tleaf.cond.fit <- fit_schoolfield.tleaf.cond_gam(data, "tleaf", "photo", "cond", T_ref = 25, start_params = start_params)

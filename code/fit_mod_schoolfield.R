schoolfield <- function(temp, J_ref, E, E_D, T_opt, T_ref = 25) {
  # temp   : Temperature values to evaluate at (C)
  # J_ref  : Rate at T_ref (units depend on rate)
  # E      : Activation energy (eV; E > 0)
  # E_D    : Deactivation energy (eV; E_D > 0)
  # T_opt  : Optimum or peak T (C)
  # T_ref  : Reference temperature for normalization (C)
  k <- 8.62e-5            # Boltzmann's constant (eV K-1)
  temp = temp + 273.15    # Convert to Kelvin
  T_ref = T_ref + 273.15  # Convert to Kelvin
  T_opt = T_opt + 273.15  # Convert to Kelvin
  # Evaluate and return
  return( J_ref * exp( E * (1/(k*T_ref) - 1/(k*temp)) ) / ( 1 + E/(E_D - E) * exp( (E_D/k) * (1/T_opt - 1/temp) ) ) )
}

fit_mod_schoolfield <- function(Data, curveID = "curveID", x = "Tleaf", y = "A", T_ref = 25, fitted_models = TRUE) {
  models <- c("schoolfield")
  aic_values <- list()
  fits <- list()  # Store model fits
  Data$T_ref = T_ref
  
  # Create a PDF file to save the plots
  pdf("schoolfield_fits.pdf")
  
  # Standardize names of x and y for use throughout
  if (x != "Tleaf") {
    Data$Tleaf = Data[[x]] # Rename X as Tleaf
    Data[[x]] = NULL      # Get rid of originals
  }
  if (y != "A") {
    Data$A = Data[[y]] # Rename Y as A
    Data[[y]] = NULL      # Get rid of originals
  }
  
  Data = subset(Data, A > 0)
  
  #Data = assign_curve_type(Data)
  Data$curve_type = "Schoolfield"
  Data = Data %>% arrange(by_group = curveID)
  results = c()
  
  # Check for missing values
  if (anyNA(Data$Tleaf) || anyNA(Data$A)) {
    stop("Missing values detected in the input data.")
  }
  # Check lengths
  if (length(Data$Tleaf) != length(Data$A)) {
    stop("The lengths of Tleaf and A columns do not match.")
  }
  # Initialize a list to store the results
  schoolfield.results.data <- list()
  
  ######################
  # Curve fitting loop #
  ######################
  # Curve fitting loop
  for (curve_id in unique(Data[[curveID]])) {
    curve_data <- subset(Data, curveID == curve_id) # Subset data for the current curveID
    print(paste("Curve ID:", curve_id, ", Number of observations:", nrow(curve_data)))  # Debug
    
    # Create a new page in the PDF for each curve
    plot_counter = 0
    
    plot <- ggplot(curve_data, aes(x = Tleaf, y = A)) +
      geom_point() +
      theme_classic() +
      ggtitle(paste("Curve ID:", curve_id))
    
    # Fit Schoolfield function using nls_multstart
    for (model in models) {
      # Fit the models
      if (model == "schoolfield") {
        # Fit schoolfield
        fit <- tryCatch({
          nls.multstart::nls_multstart(A ~ schoolfield(temp = Tleaf, J_ref, E, E_D, T_opt, T_ref),
                                       data = curve_data,
                                       iter = c(2,2,2,2),
                                       start_lower = c(J_ref = 0, E = 0, E_D = 0.2, T_opt = 0),
                                       start_upper = c(J_ref = 20, E = 2, E_D = 5, T_opt = 50),
                                       supp_errors = 'Y',
                                       na.action = na.omit,
                                       lower = c(J_ref = -10, E = 0, E_D = 0, T_opt = 0))
        }, error = function(e) {
          print(paste("Fitting failed for Curve ID:", curve_id))
          print(e)
          return(NULL)
        })
        
        # Store predictions for schoolfield
        predicted_schoolfield <- fitted(fit)
        
        # Calculate R-squared for schoolfield
        if (!is.null(fit)) {
          r_sq = 1 - sum(resid(fit)^2) / sum((curve_data$A - mean(curve_data$A))^2)
        } else {
          r_sq <- NA
        }
      }
      
      if (!is.null(fit)) { # If the fit was successful, extract parameters estimates and SE
        J_ref <- coef(fit)["J_ref"]
        J_ref_SE <- summary(fit)$coefficients[,'Std. Error']['J_ref']
        E <- coef(fit)["E"]
        E_SE <- summary(fit)$coefficients[,'Std. Error']['E']
        E_D <- coef(fit)["E_D"]
        E_D_SE <- summary(fit)$coefficients[,'Std. Error']['E_D']
        T_opt <- coef(fit)["T_opt"]
        T_opt_SE <- summary(fit)$coefficients[,'Std. Error']['T_opt']
        AIC <- AIC(fit)
        
        # Predict A values for the entire range of Tleaf
        predicted_A <- predict(fit, newdata = data.frame(Tleaf = curve_data$Tleaf))
        
        # Calculate breadth based on the 50th percentile of the predicted A values
        quant_fitted_A <- quantile(predicted_A, probs = 0.5)
        # Find Tleaf values corresponding to A at 50% of max fitted A
        idx <- which(predicted_A >= quant_fitted_A)
        if (length(idx) > 1) {
          # Take the first and last indices where A >= 50th percentile
          idx_first <- min(idx)
          idx_last <- max(idx)
          # Check if both indices are within the range of observed Tleaf values
          if (idx_first > 1 && idx_last < nrow(curve_data)) {
            breadth <- curve_data$Tleaf[idx_last] - curve_data$Tleaf[idx_first]
          } else {
            breadth <- NA
            # Add annotation if extrapolation occurred
            plot <- plot + annotate("text", x = max(curve_data$Tleaf), y = max(curve_data$A), 
                                    label = "Extrapolated Breadth", hjust = 1, vjust = -2.5, color = "red")
          }
        } else {
          # Set breadth to NA if fewer than 2 points above 50th percentile
          breadth <- NA
          # Add annotation if extrapolation occurred
          plot <- plot + annotate("text", x = max(curve_data$Tleaf), y = max(curve_data$A), 
                                  label = "Insufficient Data for Breadth", hjust = 1, vjust = -2.5, color = "red")
        }
        schoolfield.results.data[[curve_id]] <- data.frame(curveID = curve_id, J_ref, J_ref_SE, E, E_SE, E_D, E_D_SE, T_opt, T_opt_SE, AIC, r_sq, breadth)
      }
      
      if (!is.null(fit)) {
        fits[[model]] <- fit
        
        # Add fitted curves to the plots
        if (!is.null(predicted_schoolfield)) {
          plot <- plot + geom_line(aes(y = predicted_schoolfield), color = "blue") +
            geom_text(aes(x = max(curve_data$Tleaf), y = max(curve_data$A), 
                          label = paste("Schoolfield R^2 =", round(r_sq, 3))), hjust = 1, vjust = -1, color = "blue")+
            geom_hline(yintercept = quant_fitted_A, linetype = "dashed", color = "red")
        }
      } else {
        # Use provided fitted models
        if (!is.null(fitted_models[[model]])) {
          fit <- fitted_models[[model]]
          # Calculate AIC values
          aic_values[[model]] <- AIC(fit)
        }
      }
    }
    print(plot)
  }
  
  
  # Combine results for all curve IDs into a single data frame
  schoolfield.results.data <- do.call(rbind, schoolfield.results.data)
  
  # Close the PDF device
  dev.off()
  
  # Return the combined data frame
  return(schoolfield.results.data)
}


#schoolfield.results.data=fit_mod_schoolfield(at.subset3%>% select(Tleaf, A, curveID), x = "Tleaf", y = "A", T_ref = 25)
#sum(schoolfield.results.data$AIC) # = -37008.26

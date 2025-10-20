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

fit_mod_schoolfield_eWUE <- function(Data, curveID = "curveID", x = "Tleaf", y = "eWUE", T_ref = 25, fitted_models = TRUE) {
  models <- c("schoolfield")
  aic_values <- list()
  fits <- list()
  Data$T_ref = T_ref
  pdf("schoolfield_eWUE_fits.pdf")
  
  if (x != "Tleaf") {
    Data$Tleaf = Data[[x]]
    Data[[x]] = NULL
  }
  if (y != "eWUE") {
    Data$eWUE = Data[[y]]
    Data[[y]] = NULL
  }
  
  Data = subset(Data, eWUE > 0)
  Data$curve_type = "Schoolfield"
  Data = Data %>% arrange(by_group = curveID)
  schoolfield.results.data <- list()
  
  for (curve_id in unique(Data[[curveID]])) {
    curve_data <- subset(Data, curveID == curve_id)
    print(paste("Curve ID:", curve_id, ", Number of observations:", nrow(curve_data)))
    
    plot <- ggplot(curve_data, aes(x = Tleaf, y = eWUE)) +
      geom_point() +
      theme_classic() +
      ggtitle(paste("Curve ID:", curve_id))
    
    for (model in models) {
      fit <- tryCatch({
        nls.multstart::nls_multstart(eWUE ~ schoolfield(temp = Tleaf, J_ref, E, E_D, T_opt, T_ref),
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
      
      if (!is.null(fit)) {
        predicted_schoolfield <- fitted(fit)
        r_sq <- 1 - sum(resid(fit)^2) / sum((curve_data$eWUE - mean(curve_data$eWUE))^2)
        
        J_ref <- coef(fit)["J_ref"]
        J_ref_SE <- summary(fit)$coefficients[,'Std. Error']['J_ref']
        E <- coef(fit)["E"]
        E_SE <- summary(fit)$coefficients[,'Std. Error']['E']
        E_D <- coef(fit)["E_D"]
        E_D_SE <- summary(fit)$coefficients[,'Std. Error']['E_D']
        T_opt <- coef(fit)["T_opt"]
        T_opt_SE <- summary(fit)$coefficients[,'Std. Error']['T_opt']
        AIC <- AIC(fit)
        
        # Predict across observed Tleaf range
        predicted_eWUE <- predict(fit, newdata = data.frame(Tleaf = curve_data$Tleaf))
        df <- data.frame(Tleaf = curve_data$Tleaf, eWUE = predicted_eWUE)
        df <- df[order(df$Tleaf), ]
        
        # --- New breadth calc: 95% of eWUEmax ---
        eWUEmax <- max(df$eWUE, na.rm = TRUE)
        eWUE95 <- 0.95 * eWUEmax
        idx_max <- which.max(df$eWUE)
        
        left_df <- df[1:idx_max, ]
        right_df <- df[idx_max:nrow(df), ]
        
        # Find T at which eWUE crosses 95% of eWUEmax on each side
        Tlower <- left_df$Tleaf[which.min(abs(left_df$eWUE - eWUE95))]
        Tupper <- right_df$Tleaf[which.min(abs(right_df$eWUE - eWUE95))]
        
        if (Tlower == min(df$Tleaf) || Tupper == max(df$Tleaf)) {
          breadth_95 <- NA
        } else {
          breadth_95 <- Tupper - Tlower
        }
        
        schoolfield.results.data[[curve_id]] <- data.frame(
          curveID = curve_id, J_ref, J_ref_SE, E, E_SE, E_D, E_D_SE,
          T_opt, T_opt_SE, AIC, r_sq, breadth_95
        )
        
        # Add annotations to plot
        plot <- plot +
          geom_line(aes(y = predicted_schoolfield), color = "blue") +
          geom_hline(yintercept = eWUE95, linetype = "dashed", color = "purple") +
          geom_vline(xintercept = c(Tlower, Tupper), linetype = "dotted", color = "purple") +
          geom_vline(xintercept = T_opt, linetype = "dashed", color = "blue") +
          annotate("text", x = max(curve_data$Tleaf), y = eWUE95,
                   label = "95% breadth", hjust = 1, vjust = -1.5, color = "purple")
      }
    }
    print(plot)
  }
  
  schoolfield.results.data <- do.call(rbind, schoolfield.results.data)
  dev.off()
  return(schoolfield.results.data)
}

#schoolfield.results.data=fit_mod_schoolfield(at.subset3%>% select(Tleaf, eWUE, curveID), x = "Tleaf", y = "eWUE", T_ref = 25)
#sum(schoolfield.results.data$eWUEIC) # = -37008.26

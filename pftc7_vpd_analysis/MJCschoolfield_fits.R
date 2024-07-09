library(ggplot2)
library(dplyr)
library(nls.multstart)

# Define the function for MJC schoolfield model
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

# Define the fitting function
fit_mod_MJCschoolfield <- function(Data, curveID = "curveID", x = "Tleaf", y = "A", T_ref = 25, Tair = "Tair", fitted_models = NULL) {
  models <- c("mjcschoolfield")
  aic_values <- list()
  fits <- list()  # Store model fits
  Data$T_ref <- T_ref
  Data$Tair <- Data[[Tair]]  # Assign Tair from Data
  
  # Create a PDF file to save the plots
  pdf("mjcschoolfield_fits.pdf")
  
  # Standardize names of x and y for use throughout
  if (x != "Tleaf") {
    Data$Tleaf <- Data[[x]] # Rename X as Tleaf
    Data[[x]] <- NULL      # Get rid of originals
  }
  if (y != "A") {
    Data$A <- Data[[y]] # Rename Y as A
    Data[[y]] <- NULL      # Get rid of originals
  }
  
  Data <- subset(Data, A > 0)
  
  # Data = assign_curve_type(Data)
  Data$curve_type <- "mjcschoolfield"
  Data <- Data %>% arrange(.by_group = FALSE, curveID)
  
  # Check for missing values
  if (anyNA(Data$Tleaf) || anyNA(Data$A)) {
    stop("Missing values detected in the input data.")
  }
  # Check lengths
  if (length(Data$Tleaf) != length(Data$A)) {
    stop("The lengths of Tleaf and A columns do not match.")
  }
  # Initialize a list to store the results
  mjcschoolfield.results.data <- list()
  
  ######################
  # Curve fitting loop #
  ######################
  for (curve_id in unique(Data[[curveID]])) {
    curve_data <- subset(Data, curveID == curve_id)
    print(paste("Curve ID:", curve_id, ", Number of observations:", nrow(curve_data)))  # Debug
    
    # Create a new page in the PDF for each curve
    plot_counter <- 0
    
    plot <- ggplot(curve_data, aes(x = Tleaf, y = A)) +
      geom_point() +
      theme_classic() +
      ggtitle(paste("Curve ID:", curve_id))
    
    # Fit mjcschoolfield function using nls_multstart
    for (model in models) {
      # Fit the models
      if (model == "mjcschoolfield") {
        # Fit mjcschoolfield
        fit <- tryCatch({
          nls.multstart::nls_multstart(A ~ mjcschoolfield(temp = Tleaf, J_ref, E, E_D, T_opt, T_ref, Tair),
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
        
        # Store predictions for mjcschoolfield
        predicted_mjcschoolfield <- if (!is.null(fit)) fitted(fit) else NULL
        
        # Calculate R-squared for mjcschoolfield
        if (!is.null(fit)) {
          r_sq <- 1 - sum(resid(fit)^2) / sum((curve_data$A - mean(curve_data$A))^2)
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
        print(paste("Curve ID:", curve_id, "J_ref:", J_ref, "E:", E, "E_D:", E_D, "T_opt:", T_opt, "AIC:", AIC, "R_sq:", r_sq, "Breadth:", breadth))  # Debug
        mjcschoolfield.results.data[[curve_id]] <- data.frame(curveID = curve_id, J_ref, J_ref_SE, E, E_SE, E_D, E_D_SE, T_opt, T_opt_SE, AIC, r_sq, breadth)
      } else {
        print(paste("Fit was NULL for Curve ID:", curve_id))
      }
      
      if (!is.null(fit)) {
        fits[[model]] <- fit
        
        # Add fitted curves to the plots
        if (!is.null(predicted_mjcschoolfield)) {
          plot <- plot + geom_line(aes(y = predicted_mjcschoolfield), color = "blue") +
            geom_text(aes(x = max(curve_data$Tleaf), y = max(curve_data$A), 
                          label = paste("mjcschoolfield R^2 =", round(r_sq, 3))), hjust = 1, vjust = -1, color = "blue") +
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
  if (length(mjcschoolfield.results.data) > 0) {
    mjcschoolfield.results.data <- do.call(rbind, mjcschoolfield.results.data)
  } else {
    mjcschoolfield.results.data <- data.frame()
  }
  
  # Close the PDF device
  dev.off()
  
  # Return the combined data frame
  return(mjcschoolfield.results.data)
}

# Call the fit_mod_MJCschoolfield function
mjc.schoolfield.test <- fit_mod_MJCschoolfield(schld.data %>% select(Tleaf, A, curveID, Tair), x = "Tleaf", y = "A", T_ref = 25, Tair = "Tair")

write.csv(mjc.schoolfield.test, "MJCschoolfield.fits.csv")








schld.data = read.csv("outputs/raw.discardHooks_data.csv")
schld.test = schld.data%>%
  filter(curveID==c(3,6,1162))

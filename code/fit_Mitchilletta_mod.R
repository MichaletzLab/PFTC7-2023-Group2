#Load libraries
library(nls.multstart)
library(broom)
library(dplyr)
library(rTPC)
library(ggplot2)

Mitchilletta_2009 <- function (temp, b, topt, omega) {
  inside <- ((temp-topt)/b)*pi
  #omega = thermal breadth
  r_t = (1/(2*b))*(1+cos(inside))*omega
  return(r_t)
}

fit_mitchilletta <- function(Data, curveID = "curveID", x = "Tleaf", y = "A", fitted_models = NULL) {
  models <- c("Mitchilletta_2009")
  aic_values <- list()
  fits <- list()  # Store model fits
  
  # Create a PDF file to save the plots
  pdf("mitchilletta_fits.pdf")
  
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
  # Check for missing values
  if (anyNA(Data$Tleaf) || anyNA(Data$A)) {
    stop("Missing values detected in the input data.")
  }
  # Check lengths
  if (length(Data$Tleaf) != length(Data$A)) {
    stop("The lengths of Tleaf and A columns do not match.")
  }
  # Initialize a list to store the results
  mitchilletta.results.data <- list()
  
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
    
    for (model in models) {
      # Fit the models
      if (model == "Mitchilletta_2009") {
        # Fit weibull_1995
        fit <- tryCatch({
          nls.multstart::nls_multstart(A ~ Mitchilletta_2009(temp = Tleaf, b, topt, omega),
                                       data = curve_data,
                                       iter = c(3,3,3),
                                       start_lower = c(b = 0, topt = 10, omega = 1),
                                       start_upper = c(b = 50, topt = 50, omega = 50),
                                       lower = c(b=0, topt=0, omega=0),
                                       upper = c(b=50, topt=50, omega=90),
                                       supp_errors = 'Y',
                                       convergence_count = FALSE)
        }, error = function(e) {
          print("BadFit")
          return(NULL)
        })
        # Store predictions for weibull_1995
        predicted_mitchilletta <- fitted(fit)
        # Calculate R-squared for weibull_1995
        RSS <- sum(residuals(fit)^2) # Calculate residual sum of squares
        mean_y <- mean(curve_data$A) # Calculate total sum of squares
        TSS <- sum((curve_data$A - mean_y)^2)
      }
      
      if(typeof(fit) != "NULL") { # If the fit was successful, extract parameters estimates and SE
        b <- coef(fit)["b"]
        b_SE <- summary(fit)$coefficients[,'Std. Error']['b']
        T_opt <- coef(fit)["topt"]
        T_opt_SE <- summary(fit)$coefficients[,'Std. Error']['topt']
        omega <- coef(fit)["omega"]
        omega_SE <- summary(fit)$coefficients[,'Std. Error']['omega']
        AIC <- AIC(fit)
        r_sq <- 1 - sum(resid(fit)^2) / sum((curve_data$A - mean(curve_data$A))^2)
        
        # Store the results for this curve ID
        mitchilletta.results.data[[curve_id]] <- data.frame(curveID=curve_id, b, b_SE, T_opt, T_opt_SE,omega, omega_SE, AIC, r_sq)
      }
      if (!is.null(fit)) {
        fits[[model]] <- fit
        
        # Add fitted curves to the plots
        if (!is.null(predicted_mitchilletta)) {
          plot <- plot + geom_line(aes(y = predicted_mitchilletta), color = "blue") +
            geom_text(aes(x = max(curve_data$Tleaf), y = max(curve_data$A), 
                          label = paste("mitchilletta R^2 =", round(r_sq, 3))), hjust = 1, vjust = -1, color = "blue")
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
  mitchilletta.results.data <- do.call(rbind, mitchilletta.results.data)
  
  # Close the PDF device
  dev.off()
  
  # Return the combined data frame
  return(mitchilletta.results.data)
}


mitchilletta.results.data <- fit_mitchilletta(at.subset3)
sum(mitchilletta.results.data$AIC) # = 141976.7
#write.csv(weibull.results.data, file.path("data", "weibull_results.csv"), row.names = FALSE)

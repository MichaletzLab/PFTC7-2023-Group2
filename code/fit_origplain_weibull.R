#Load libraries
library(nls.multstart)
library(broom)
library(dplyr)
library(rTPC)
library(ggplot2)

#The weibull90 from padfield calculates b as the breadth at 50% of Amax. I figured this out manually.
#This function modifies the original to calculate breadth as 90% of Amax. 
#This can be changed by modifying the ln expression located after the b in the equation. 
# Since the original is Q(50) we have to do Q(90) / Q(50)
weibull90 <- function (temp, a, topt, b, c) {
  k <- 8.62e-5
  cexpn1 <- (((c-1)/c)^((1-c)/c))
  cexpn2 <- (((c-1)/c)^(1/c))
  cexpn3 <- ((c-1)/c)
  #newfrac <- (((-log(0.1))^(1/k)) / ((-log(0.5))^(1/k))) #This is Q(90) / Q(50)
  newfrac <- log(0.1)/log(0.5)
  tempexpn <- ((temp-topt)/(b*newfrac))
  oldexpn <- ((temp-topt)/(b))
  rt <- a * cexpn1 * ((oldexpn+cexpn2)^(c-1)) * (exp(-(oldexpn+cexpn2)^c) + cexpn3)
  return(rt)
}

fit_weibull <- function(Data, curveID = "curveID", x = "Tleaf", y = "A") {
  models <- c("weibull90")
  aic_values <- list()
  fits <- list()  # Store model fits
  
  # Create a PDF file to save the plots
  pdf("weibull_fits.pdf")
  
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
  weibull90.results.data <- list()
  
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
      if (model == "weibull90") {
        # Fit weibull90
        fit <- tryCatch({
          nls.multstart::nls_multstart(A ~ weibull90(temp = Tleaf, a, topt, b, c),
                                       data = curve_data,
                                       iter = c(3,3,3,3),
                                       start_lower = c(a = 0,topt = 10, b = 0, c = 0),
                                       start_upper = c(a = 50, topt = 50, b = 30, c = 30),
                                       lower = c(a = 0,topt = 0, b = 0, c = 0),
                                       upper = c(a = 50, topt = 50, b = 50, c = 50),
                                       supp_errors = 'Y',
                                       convergence_count = FALSE)
        }, error = function(e) {
          print("BadFit")
          return(NULL)
        })
        # Store predictions for weibull90
        predicted_weibull <- fitted(fit)
        # Calculate R-squared for weibull90
        RSS <- sum(residuals(fit)^2) # Calculate residual sum of squares
        mean_y <- mean(curve_data$A) # Calculate total sum of squares
        TSS <- sum((curve_data$A - mean_y)^2)
      }
      
      if(typeof(fit) != "NULL") { # If the fit was successful, extract parameters estimates and SE
        a <- coef(fit)["a"]
        a_SE <- summary(fit)$coefficients[,'Std. Error']['a']
        b <- coef(fit)["b"]
        b_SE <- summary(fit)$coefficients[,'Std. Error']['b']
        c_D <- coef(fit)["c"]
        c_SE <- summary(fit)$coefficients[,'Std. Error']['c']
        T_opt <- coef(fit)["topt"]
        T_opt_SE <- summary(fit)$coefficients[,'Std. Error']['topt']
        AIC <- AIC(fit)
        r_sq <- 1 - sum(resid(fit)^2) / sum((curve_data$A - mean(curve_data$A))^2)
        
        # Store the results for this curve ID
        weibull.results.data[[curve_id]] <- data.frame(curveID=curve_id, a, a_SE, b, b_SE, c_D, c_SE, T_opt, T_opt_SE, AIC, r_sq)
      }
      if (!is.null(fit)) {
        fits[[model]] <- fit
        
        # Add fitted curves to the plots
        if (!is.null(predicted_weibull)) {
          plot <- plot + geom_line(aes(y = predicted_weibull), color = "blue") +
            geom_text(aes(x = max(curve_data$Tleaf), y = max(curve_data$A), 
                          label = paste("Weibull R^2 =", round(r_sq, 3))), hjust = 1, vjust = -1, color = "blue")
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
  weibull.results.data <- do.call(rbind, weibull.results.data)
  
  # Close the PDF device
  dev.off()
  
  # Return the combined data frame
  return(weibull.results.data)
}

weibull90.results.data <- fit_weibull(at.subset3)

#sum(weibull90.results.data$AIC) # = -32437.7
#weibull90.results.data <- fit_weibull(at.subset3 %>% select(Tleaf, A, curveID), x = "Tleaf", y = "A")
#write.csv(weibull90.results.data, file.path("data", "weibull_results.csv"), row.names = FALSE)

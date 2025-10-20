#Load libraries
library(nls.multstart)
library(broom)
library(dplyr)
library(rTPC)
library(ggplot2)

fit_weibull.gsw_breadth <- function(Data, curveID = "curveID", x = "Tleaf", y = "gsw") {
  models <- c("weibull_1995")
  aic_values <- list()
  fits <- list()  # Store model fits
  
  # Create a PDF file to save the plots
  pdf("weibull.gsw_fits.pdf")
  
  # Standardize names of x and y for use throughout
  if (x != "Tleaf") {
    Data$Tleaf <- Data[[x]] # Rename X as Tleaf
    Data[[x]] <- NULL      # Get rid of originals
  }
  if (y != "gsw") {
    Data$gsw <- Data[[y]] # Rename Y as A
    Data[[y]] <- NULL      # Get rid of originals
  }
  Data <- subset(Data, gsw > 0)
  # Check for missing values
  if (anyNA(Data$Tleaf) || anyNA(Data$gsw)) {
    stop("Missing values detected in the input data.")
  }
  # Check lengths
  if (length(Data$Tleaf) != length(Data$gsw)) {
    stop("The lengths of Tleaf and A columns do not match.")
  }
  # Initialize a list to store the results
  weibull.gsw.results.data <- list()
  prediction.data <- list()
  
  # Curve fitting loop
  for (curve_id in unique(Data[[curveID]])) {
    curve_data <- subset(Data, curveID == curve_id) # Subset data for the current curveID
    print(paste("Curve ID:", curve_id, ", Number of observations:", nrow(curve_data)))  # Debug
    
    # Create a new page in the PDF for each curve
    plot_counter <- 0
    
    plot <- ggplot(curve_data, aes(x = Tleaf, y = gsw)) +
      geom_point() +
      theme_classic() +
      ggtitle(paste("Curve ID:", curve_id))
    
    for (model in models) {
      # Fit the models
      if (model == "weibull_1995") {
        # Fit weibull_1995
        fit <- tryCatch({
          nls.multstart::nls_multstart(gsw ~ weibull_1995(temp = Tleaf, a, topt, b, c),
                                       data = curve_data,
                                       iter = c(3,3,3,3),
                                       start_lower = c(a = 0,topt = 10, b = 0, c = 0),
                                       start_upper = c(a = 50, topt = 50, b = 30, c = 30),
                                       lower = get_lower_lims(curve_data$Tleaf, curve_data$A, model_name = 'weibull_1995'),
                                       upper = get_upper_lims(curve_data$Tleaf, curve_data$A, model_name = 'weibull_1995'),
                                       supp_errors = 'Y',
                                       convergence_count = FALSE)
        }, error = function(e) {
          print("BadFit")
          return(NULL)
        })
        # Store predictions for weibull_1995
        predicted_weibull <- fitted(fit)
        # Calculate R-squared for weibull_1995
        RSS <- sum(residuals(fit)^2) # Calculate residual sum of squares
        mean_y <- mean(curve_data$gsw) # Calculate total sum of squares
        TSS <- sum((curve_data$gsw - mean_y)^2)
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
        
        predicted_gsw <- predict(fit, newdata = data.frame(Tleaf = curve_data$Tleaf))
        
        # Store the predicted values along with Tleaf and curveID
        prediction.data[[curve_id]] <- data.frame(curveID = curve_id, Tleaf = curve_data$Tleaf, predicted_gsw)
        
        # Calculate Amax
        gswmax <- max(predicted_gsw)
        # Calculate breadth at 80%
        getbreadth_90 <- get_breadth(fit, level = 0.95)
        
        # Try to get breadth another way
        gswmin <- min(predicted_gsw)
        gswrange <- gswmax - gswmin
        gsw90max <- 0.95 * gswrange
        gswat90b <- gswmin + gsw90max
        
        # Find the two Tleaf values closest to Aat90b on either side of Amax
        df <- data.frame(Tleaf = curve_data$Tleaf, gsw = predicted_gsw)
        df <- df[order(df$Tleaf), ]
        idx_max <- which.max(df$gsw)
        
        # Find closest value to Aat90b on the left of Amax
        left_df <- df[1:idx_max, ]
        Tlower <- left_df$Tleaf[which.min(abs(left_df$gsw - gswat90b))]
        
        # Find closest value to Aat90b on the right of Amax
        right_df <- df[idx_max:length(df$gsw), ]
        Tupper <- right_df$Tleaf[which.min(abs(right_df$gsw - gswat90b))]
        
        # Check if Tlower and Tupper are at the data boundaries
        if (Tlower == min(df$Tleaf) || Tupper == max(df$Tleaf)) {
          mjcbreadth <- NA
        } else {
          mjcbreadth <- Tupper - Tlower
        }
        
        get_90_gsw <- gswmax*0.95
        # Check if the A value of the purple dot is below either side of the range
        if (get_90_gsw < min(df$gsw[left_df$gsw <= get_90_gsw]) || get_90_gsw < min(df$gsw[right_df$gsw >= get_90_gsw])) {
          getbreadth_90 <- NA
        }
        
        # Store the results for this curve ID
        weibull.gsw.results.data[[curve_id]] <- data.frame(curveID = curve_id, a, a_SE, b, b_SE, c_D, c_SE, T_opt, T_opt_SE, AIC, r_sq, Tlower, Tupper, mjcbreadth, getbreadth_90)
      }
      if (!is.null(fit)) {
        fits[[model]] <- fit
        
        # Add fitted curves to the plots
        if (!is.null(predicted_weibull)) {
          plot <- plot + geom_line(aes(y = predicted_weibull), color = "blue") +
            geom_text(aes(x = max(curve_data$Tleaf), y = max(curve_data$gsw), 
                          label = paste("Weibull R^2 =", round(r_sq, 3))), hjust = 1, vjust = -1, color = "blue") +
            geom_point(data = data.frame(Tleaf = c(Tlower, Tupper), gsw = predict(fit, newdata = data.frame(Tleaf = c(Tlower, Tupper)))),
                       aes(x = Tleaf, y = gsw), color = "pink", size = 3) +
            geom_point(data = data.frame(Tleaf = T_opt, gsw = get_90_gsw),
                       aes(x = Tleaf, y = gsw), color = "purple", size = 3)
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
  weibull.gsw.results.data <- do.call(rbind, weibull.gsw.results.data)
  
  # Combine prediction data for all curve IDs into a single data frame
  prediction.data <- do.call(rbind, prediction.data)
  
  # Save the prediction data to a CSV file
  write.csv(prediction.data, "prediction.gsw_data.csv", row.names = FALSE)
  
  # Close the PDF device
  dev.off()
  
  # Return the combined data frame
  return(weibull.gsw.results.data)
}

#weibull.results.data <- fit_weibull_breadth(at.subset3)

#sum(weibull_1995.results.data$AIC) # = -32437.7
#weibull_1995.results.data <- fit_weibull(at.subset3 %>% select(Tleaf, A, curveID), x = "Tleaf", y = "A")
#write.csv(weibull_1995.results.data, file.path("data", "weibull_results.csv"), row.names = FALSE)


###note you have to run this twice. Once to get the main pdf with one plot per page
  #and a second time to get the plot summary document. The info to change has been noted.

#Load libraries
library(nls.multstart)
library(broom)
library(dplyr)
library(rTPC)
library(ggplot2)
library(segmented)
library(cowplot)

schoolfield <- function(temp, J_ref, EA, E_D, T_opt, T_ref = 25) {
  # temp   : Temperature values to evaluate at (C)
  # J_ref  : Rate at T_ref (units depend on rate)
  # EA      : Activation energy (eV; E > 0)
  # E_D    : Deactivation energy (eV; E_D > 0)
  # T_opt  : Optimum or peak T (C)
  # T_ref  : Reference temperature for normalization (C)
  k <- 8.62e-5            # Boltzmann's constant (eV K-1)
  temp = temp + 273.15    # Convert to Kelvin
  T_ref = T_ref + 273.15  # Convert to Kelvin
  T_opt = T_opt + 273.15  # Convert to Kelvin
  # Evaluate and return
  return( J_ref * exp( EA * (1/(k*T_ref) - 1/(k*temp)) ) / ( 1 + EA/(E_D - EA) * exp( (E_D/k) * (1/T_opt - 1/temp) ) ) )
}
#####################################
# Curve fitting comparison function
#####################################

fit_schoolfield_weibull <- function(Data, curveID = "curveID", x = "Tleaf", y = "A", T_ref = 25, fitted_models = NULL) {
  models <- c("schoolfield", "weibull_1995")
  fits <- list()  # Store model fits
  Data$T_ref <- T_ref
  
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
  Data$curve_type <- "Schoolfield"
  Data <- Data %>% arrange(curveID)  # Change "by_group" to "curveID"
  
  # Check for missing values
  if (anyNA(Data$Tleaf) || anyNA(Data$A)) {
    stop("Missing values detected in the input data.")
  }
  
  # Create a PDF file to save the plots
  pdf("schoolfield_weibull_fits.pdf")
  plot_list <- list()
  
  # Curve fitting loop
  for (curve_id in unique(Data[[curveID]])) {
    curve_data <- subset(Data, curveID == curve_id) # Subset data for the current curveID
    
    # Initialize plot
    plot <- ggplot(curve_data, aes(x = Tleaf, y = A)) +
      geom_point() +
      theme_classic() +
      ggtitle(paste("Curve ID:", curve_id)) +
      theme(plot.margin = margin(1, 1, 1, 1, "cm")) # Add margin
    
    schoolfield_r_sq <- NA
    schoolfield_aic <- NA
    weibull_r_sq <- NA
    weibull_aic <- NA
    predicted_schoolfield <- NULL
    predicted_weibull <- NULL
    
    # Fit Schoolfield and Weibull functions using nls_multstart
    for (model in models) {
      if (model == "schoolfield") {
        # Fit schoolfield
        fit <- tryCatch({
          nls.multstart::nls_multstart(A ~ schoolfield(temp = Tleaf, J_ref, EA, E_D, T_opt, T_ref),
                                       data = curve_data,
                                       iter = c(2,2,2,2),
                                       start_lower = c(J_ref = 0, EA = 0, E_D = 0.2, T_opt = 0),
                                       start_upper = c(J_ref = 20, EA = 2, E_D = 5, T_opt = 50),
                                       supp_errors = 'Y',
                                       na.action = na.omit,
                                       lower = c(J_ref = -10, EA = 0, E_D = 0, T_opt = 0))
        }, error = function(e) {
          print(paste("Fitting failed for Curve ID:", curve_id))
          print(e)
          return(NULL)
        })
        
        if (!is.null(fit)) {
          # Predict Schoolfield values
          predicted_schoolfield <- data.frame(Tleaf = curve_data$Tleaf, A = predict(fit, newdata = curve_data))
          schoolfield_r_sq <- 1 - sum(resid(fit)^2) / sum((curve_data$A - mean(curve_data$A))^2)
          schoolfield_aic <- round(AIC(fit), 3)
        }
      } else if (model == "weibull_1995") {
        # Fit Weibull
        fit <- tryCatch({
          nls.multstart::nls_multstart(A ~ weibull_1995(temp = Tleaf, a, topt, b, c),
                                       data = curve_data,
                                       iter = c(3,3,3,3),
                                       start_lower = c(a = 0, topt = 10, b = 0, c = 0),
                                       start_upper = c(a = 50, topt = 50, b = 30, c = 30),
                                       lower = get_lower_lims(curve_data$Tleaf, curve_data$A, model_name = 'weibull_1995'),
                                       upper = get_upper_lims(curve_data$Tleaf, curve_data$A, model_name = 'weibull_1995'),
                                       supp_errors = 'Y',
                                       convergence_count = FALSE)
        }, error = function(e) {
          print("BadFit")
          return(NULL)
        })
        
        if (!is.null(fit)) {
          # Predict Weibull values
          predicted_weibull <- data.frame(Tleaf = curve_data$Tleaf, A = predict(fit, newdata = curve_data))
          weibull_r_sq <- 1 - sum(resid(fit)^2) / sum((curve_data$A - mean(curve_data$A))^2)
          weibull_aic <- round(AIC(fit), 3)
        }
      }
    }
    ###########Un-comment out the annotation below (two chunks)to produce the proper schoolfield_weibull doc
    # Add the predictions to the plot
    if (!is.null(predicted_schoolfield)) {
      plot <- plot + geom_line(data = predicted_schoolfield, aes(y = A, color = "Schoolfield")) +
        theme(legend.position = "none")+  # Remove legend
        scale_color_manual(values = c("red", "blue"), labels = c("Schoolfield", "Weibull"))#+
        #annotate("text", x = max(curve_data$Tleaf), y = min(curve_data$A), 
                 #label = paste("Schoolfield R^2 =", round(schoolfield_r_sq, 3), ", AIC =", schoolfield_aic), 
                 #hjust = 1, vjust = -1, color = "red") 
        
    }
    if (!is.null(predicted_weibull)) {
      plot <- plot + geom_line(data = predicted_weibull, aes(y = A, color = "Weibull")) +
        theme(legend.position = "none") + # Remove legend
        scale_color_manual(values = c("red", "blue"), labels = c("Schoolfield", "Weibull"))#+
        #annotate("text", x = max(curve_data$Tleaf), y = min(curve_data$A), 
                 #label = paste("Weibull R^2 =", round(weibull_r_sq, 3), ", AIC =", weibull_aic), 
                 #hjust = 1, vjust = 1, color = "blue") 
    }
    
    # Print the plot for this curve ID
    print(plot)
    plot_list[[length(plot_list) + 1]] <- plot + labs(color = "Model") + theme(text = element_text(size = 5))+theme(plot.margin = unit(c(0,0,0,0), "cm"))+
      annotate("text", x = max(curve_data$Tleaf), y = min(curve_data$A), 
               label = paste("S.field R2 =", round(schoolfield_r_sq, 3), ", AIC =", schoolfield_aic), 
               hjust = 1, vjust = -0.2, color = "red", size=1.7)+
      annotate("text", x = max(curve_data$Tleaf), y = min(curve_data$A), 
               label = paste("Weibull R2 =", round(weibull_r_sq, 3), ", AIC =", weibull_aic), 
               hjust = 1, vjust = 0.8, color = "blue", size=1.7)
  }
  
  # Close the PDF device for individual plots
  dev.off()
  
  # Create a PDF file for the summary page
  pdf("schoolfield_weibull_summary.pdf")
  
  # Determine the number of rows and columns for the plot grid
  num_plots <- length(plot_list)
  num_per_page <- 25  # Adjust as needed
  num_pages <- ceiling(num_plots / num_per_page)
  
  for (page in 1:num_pages) {
    # Determine the range of plots to include on this page
    start_plot <- (page - 1) * num_per_page + 1
    end_plot <- min(start_plot + num_per_page - 1, num_plots)
    plot_range <- start_plot:end_plot
    
    # Select the plots for this page
    plots_on_page <- plot_list[plot_range]
    
    # Arrange the plots into a grid and save to the PDF
    summary_plot <- cowplot::plot_grid(plotlist = plots_on_page, ncol = 5)
    print(summary_plot)
  }
  
  dev.off()
}

############################
#Testing subset:
test123 <- at.subset3 %>%
  filter(curveID %in% c(3, 6, 9, 12, 15, 18, 21))
testing123 <- fit_schoolfield_weibull(at.subset3, x = "Tleaf", y = "A", T_ref = 25)
###############################################################################

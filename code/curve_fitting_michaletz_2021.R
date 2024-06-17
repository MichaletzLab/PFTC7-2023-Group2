# Curve fitting master code for activation energy project. This file includes functions describing
# Arrhenius and Sharpe-Schoolfield temperature response functions, a function which assigns each
# curveID to be fit with either Arrhenius or SS (or neither) based on specified criteria, and a
# curve-fitting function which fits the assigned curve type to each dataset using nls.multstart,
# then plots the results if desired.

# Last modified Josef Garen 2 March 2021

library(nls.multstart)
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)
library(scales)
library(gridExtra)
library(grDevices)
#dev.off()

# Returns Sharpe-Schoolfield function for the given parameters
michaletz_2021 <- function(temp, J_ref, E, E_D, T_opt, T_ref = 25) {
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



# Returns Arrhenius function for the given parameters
arrhenius <- function(temp, J_ref, E, T_ref = 25) {
  # temp   : Temperature values to evaluate function at (C)
  # J_ref  : Rate at T_ref (units depend on rate)
  # E      : Activation energy (eV; E > 0)
  # T_ref  : Reference temperature for normalization (C)

  k <- 8.62e-5           # Boltzmann's constant (eV K-1)
  temp = temp + 273.15   # Convert to Kelvin
  T_ref = T_ref + 273.15 # Convert to Kelvin

  # Evaluate and return
  return( J_ref * exp( E * (1/(k*T_ref) - 1/(k*temp)) ) )
}



# In the given dataset, assigns each curveID to be fit with Arrhenius, Schoolfield, or neither
assign_curve_type <- function(data) {

  filter_criteria = c() # Holder

  #new_dat = c() ###
  # Find metrics for filtering for each curve
  for (i in unique(data$curveID)) {

    cur_dat = subset(data, curveID == i) # Grab the current data set
    Tpk = cur_dat$Tleaf[which( cur_dat$Photo == max(cur_dat$Photo) )[1]] # Find peak
    #cur_dat = subset(cur_dat, Tleaf <= Tpk) ###

    nbefore = sum( cur_dat$Tleaf < Tpk ) # Find number of points before peak
    nafter = sum( cur_dat$Tleaf > Tpk )  # Find number of points before peak
    Tbefore = Tpk-min(cur_dat$Tleaf)     # Find temperature span before peak

    filter_criteria = rbind(filter_criteria, data.frame(curveID = i, nbefore, nafter, Tbefore))

    #new_dat = rbind(new_dat, cur_dat) ###
  }

  # Merge filtering metrics with dataset for ease of filtering
  data = merge(data, filter_criteria, by="curveID")

  data$curve_type = NA
  data$failure_status = NA

  # Reject curves if the are decreasing only, have too little T span, or too few points
  data$failure_status = ifelse(data$nbefore < 2, "Rejected: fewer than two points before peak", data$failure_status)
  data$failure_status = ifelse(data$Tbefore < 4, "Rejected: ascending region too small", data$failure_status)
  data$failure_status = ifelse(data$nbefore+data$nafter < 4, "Rejected: fewer than 5 data points", data$failure_status)

  # Assign curve type based on whether or not they have an identifiable peak, or NA if failed above
  data$curve_type = ifelse(data$nbefore > 0 & data$nafter > 0, "Schoolfield", data$curve_type)
  data$curve_type = ifelse(data$nafter < 1, "Arrhenius", data$curve_type)
  data$curve_type = ifelse(!is.na(data$failure_status), NA, data$curve_type)

  data = data %>% select(-nbefore, -nafter, -Tbefore)  # Remove junk
  return(data)
}



# Fit curves, plot, and return dataframe with best fit parameters
fit_curves_michaletz_2021 <- function(Data, x = "Tleaf", y = "Photo", T_ref = 25, PLOT=TRUE) {

  set.seed(1)
  Data$T_ref = T_ref

  # Standardize names of x and y for use throughout
  if(x != "Tleaf") {
    Data$Tleaf = Data[[x]] # Rename X as Tleaf
    Data[[x]] = NULL      # Get rid of originals
  }
  if(y != "Photo") {
    Data$Photo = Data[[y]] # Rename Y as Photo
    Data[[y]] = NULL      # Get rid of originals
  }

  Data = subset(Data, Photo > 0)

  # Open output file if plotting desired
  if(PLOT==TRUE)
    grDevices::pdf("plots_all.pdf", width = 10, height = 5)

  Data = assign_curve_type(Data)
  Data = Data %>% arrange(by_group = curveID)
  results = c()

  ######################
  # Curve fitting loop #
  ######################
  for (i in unique(Data$curveID)) {

    print(paste(i, "/", as.character(max(Data$curveID)), "complete"))

    d_1 <- subset(Data, curveID == i)    # Grab only current curveID
#    d_1 <-subset(d_1, Photo > 0)         # Drop all values with Y < 0

    if (is.na(d_1$curve_type)) {     # If no curve should be fit, all NAs
#      lnB0 = lnB0_SE = E = E_SE = E_D = E_D_SE = T_h = T_h_SE = AIC = r_sq = NA
      J_ref = J_ref_SE = E = E_SE = E_D = E_D_SE = T_opt = T_opt_SE = AIC = r_sq = NA


    } else if(dim(d_1)[1] < 5) {     # If not enough positive points, all NAs and set failure status
#      lnB0 = lnB0_SE = E = E_SE = E_D = E_D_SE = T_h = T_h_SE = AIC = r_sq = NA
      J_ref = J_ref_SE = E = E_SE = E_D = E_D_SE = T_opt = T_opt_SE = AIC = r_sq = NA
      d_1$failure_status = "Rejected: fewer than 5 data points with Y > 0"
      d_1$curve_type = NA

    } else if(d_1$curve_type == "Schoolfield") { # If curve type is Schoolfield, fit one

      # Fit Schoolfield function using nls_multstart
      fit <- nls_multstart(Photo ~ michaletz_2021(temp = Tleaf, J_ref, E, E_D, T_opt, T_ref),
                           data = d_1,
                           iter = 1000,
                           start_lower = c(J_ref = 0, E = 0, E_D = 0.2, T_opt = 0),
                           start_upper = c(J_ref = 20, E = 2, E_D = 5, T_opt = 50),
                           supp_errors = 'Y',
                           na.action = na.omit,
                           lower = c(J_ref = -10, E = 0, E_D = 0, T_opt = 0))

      if(typeof(fit) != "NULL") { # If the fit was successful, extract parameters estimates and SE
        J_ref <- coef(fit)["J_ref"]
        J_ref_SE <- summary(fit)$coefficients[,'Std. Error']['J_ref']
        E <- coef(fit)["E"]
        E_SE <- summary(fit)$coefficients[,'Std. Error']['E']
        E_D <- coef(fit)["E_D"]
        E_D_SE <- summary(fit)$coefficients[,'Std. Error']['E_D']
        T_opt <- coef(fit)["T_opt"]
        T_opt_SE <- summary(fit)$coefficients[,'Std. Error']['T_opt']
        AIC <- AIC(fit)
        r_sq = 1-sum(resid(fit)^2)/sum((d_1$Photo-mean(d_1$Photo))^2)

        if (T_opt < min(d_1$Tleaf)) {
          J_ref = J_ref_SE = E = E_SE = E_D = E_D_SE = T_opt = T_opt_SE = AIC = r_sq = NA
          d_1$failure_status = "Poor fit: T_opt out of bounds"
        }

      } else { # Otherwise, likely a convergence failure. All set to NA
 #       lnB0 = lnB0_SE = E = E_SE = E_D = E_D_SE = T_h = T_h_SE = AIC = r_sq = NA
        J_ref = J_ref_SE = E = E_SE = E_D = E_D_SE = T_opt = T_opt_SE = AIC = r_sq = NA
        d_1$failure_status = "Convergence failure"
      }

    } else if(d_1$curve_type == "Arrhenius") { # If curve type is Arrhenius, fit one

      # Fit Arrhenius function using nls_multstart
      fit <- nls_multstart(Photo ~ arrhenius(temp = Tleaf, J_ref, E, T_ref = T_ref),
                           data = d_1,
                           iter = 1000,
                           start_lower = c(J_ref = 0, E = 0),
                           start_upper = c(J_ref = 20, E = 4),
                           supp_errors = 'Y',
                           na.action = na.omit,
                           lower = c(lnB0 = -10, E = 0))


      if(typeof(fit) != "NULL") { # If the fit was successful, extract parameters estimates and SE
        J_ref <- coef(fit)["J_ref"]
        J_ref_SE <- summary(fit)$coefficients[,'Std. Error']['J_ref']
        E <- coef(fit)["E"]
        E_SE <- summary(fit)$coefficients[,'Std. Error']['E']
        E_D = E_D_SE = T_opt = T_opt_SE = NA
        AIC <- AIC(fit)
        r_sq = 1-sum(resid(fit)^2)/sum((d_1$Photo-mean(d_1$Photo))^2)

      } else { # Otherwise, likely a convergence failure. All set to NA
#        lnB0 = lnB0_SE = E = E_SE = E_D = E_D_SE = T_h = T_h_SE = AIC = r_sq = NA
        J_ref = J_ref_SE = E = E_SE = E_D = E_D_SE = T_opt = T_opt_SE = AIC = r_sq = NA
        d_1$failure_status = "Convergence failure"
      }
    }


    ##########################
    # Plot curves if desired #
    ##########################
    if (PLOT==TRUE) {

      # Pick a set of temperature values to compute model over
      tmp_temps <- seq(min(floor(d_1$Tleaf)), ceiling(max(d_1$Tleaf)), length = 200)

      # Compute model predictions depending on model type and select curve color
      if(is.na(d_1$curve_type[1])) {
        tmp_model = mean(d_1$Photo)
        color_model = "white"
      } else if(d_1$curve_type[1] == "Schoolfield") {
        tmp_model <- michaletz_2021(tmp_temps, J_ref, E, E_D, T_opt, T_ref = T_ref)
        color_model = "red"
      } else if(d_1$curve_type[1] == "Arrhenius") {
        tmp_model <- arrhenius(tmp_temps, J_ref, E, T_ref = T_ref)
        color_model = "blue"
      }

      # Assemble model predicitons and raw data to plot
      ModelToPlotS <- data.frame(Temperature = tmp_temps, TraitValue = tmp_model)
      DataToPlot <- na.omit(data.frame(Temperature = d_1$Tleaf, TraitValue = d_1$Photo))

      # Build plots - normal
      p <- ggplot() +
        geom_point(data = DataToPlot, aes(x = Temperature, y = TraitValue),
                   shape = 21, size = 2.75, col = "black", fill = "white") +
        geom_line(data = ModelToPlotS,aes(x = Temperature, y = TraitValue), colour = color_model) +
        ggtitle(paste("curveID=", i, ", ", d_1$Taxon[1], sep="")) +
        xlab(expression(paste("Leaf temperature (", degree, C, ")"))) +
        ylab(expression("Assimilation rate (" * mu ~ "mol" ~m^-2 ~s^-1 * ")")) +
        theme_bw(base_size=12) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black")) +
        annotate("text", size = 3,
                 label=paste("R^2 = ", sprintf("%.2f", r_sq),"\nE =", format(E, digits = 3),"\nAIC =",format(AIC,digits=3), "\nCurve type = ", d_1$curve_type[1]),
                 x = min(DataToPlot[, "Temperature"]),
                 y = mean(DataToPlot[, "TraitValue"]),
                 hjust=0,
                 fontface = 3)

      # Build plots - Arrhenius plot
      p3 <- ggplot() +
#        geom_point(data = DataToPlot, aes(x = 1/(0.00008617*(Temperature + 273.15)), y = TraitValue),
#                   shape = 21, size = 2.75, col = "black", fill = "white") +
#        geom_line(data = ModelToPlotS, aes(x = 1/(0.00008617*(Temperature + 273.15)), y = TraitValue),
#                  colour = color_model) +
        geom_point(data = DataToPlot, aes(x = 1/(0.00008617*(T_ref + 273.15)) - 1/(0.00008617*(Temperature + 273.15)), y = TraitValue),
                   shape = 21, size = 2.75, col = "black", fill = "white") +
        geom_line(data = ModelToPlotS, aes(x = 1/(0.00008617*(T_ref + 273.15)) - 1/(0.00008617*(Temperature + 273.15)), y = TraitValue),
                  colour = color_model) +
        ggtitle(paste("curveID=", i, ", ", d_1$Taxon[1], sep="")) +
#        xlab(expression(paste('Leaf temperature ', '1/',italic('kT'), ' (',  eV^{-1}, ')'))) +
        xlab(expression(paste('Leaf temperature ', '[ 1/',italic('kT_ref'), ' - 1/',italic('kT'), '] (',  eV^{-1}, ')'))) +
        ylab(expression("Assimilation rate (" * mu ~ "mol" ~m^-2 ~s^-1 * ")")) +
#        scale_x_continuous(sec.axis = sec_axis(trans = ~ (1/(.*0.00008617))-273.15 , name = expression(paste("Leaf temperature (", degree, C, ")")))) +
        scale_x_continuous(sec.axis = sec_axis(trans = ~ (T_ref+273.15)/(1 - 0.00008617*(T_ref+273.15)*.) - 273.15, name = expression(paste("Leaf temperature (", degree, C, ")")))) +
        scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3),
                           labels = trans_format("log", math_format(e^.x))) +
        theme_bw(base_size=12) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"))

      grid.arrange(p,p3,nrow=1) # Arrange plots next to each other
    } # End plotting

    # Build dataframe with results
    results = rbind(results, data.frame(d_1[1,] %>% select(-Tleaf, -Photo), J_ref, J_ref_SE, E, E_SE,
                                        E_D, E_D_SE, T_opt, T_opt_SE, AIC, r_sq))
  }

  # Close output file
  if(PLOT==TRUE)
    grDevices::dev.off()

  return(results)
}


# Fit quadratic AT response


quadratic.model = function(temp, a, b, c) {

  return(a*temp*temp + b*temp + c)

}

fit.quadratic = function(data.AT) {


  output.params = c()

  for (i in unique(data.AT$curveID)) {
    print(i)
    d_1 = subset(data.AT, curveID == i)

    fit = nls_multstart(A ~ quadratic.model(temp = Tleaf, a1,b1,c1),
                        data = d_1,
                        iter = 1000,
                        start_lower =  c(a1 = -5, b1 = 0, c1 = -10),
                        start_upper = c(a1 = 0, b1 = 10, c1 = 10),
                        supp_errors = 'Y',
                        na.action = na.omit)


    #fit = lm(A ~ poly(Tleaf, degree = 2, raw = T), data = d_1)

    # If the fit failed, set failure status and all params to NA
    if(typeof(fit) == "NULL" || coef(fit)["a1"] >= 0) {
      Topt = Amax = Tbreadth = CTmin = CTmax = r_sq = NA
      print(coef(fit)["a1"])

    } else {

      a <- coef(fit)["a1"]
      b <- coef(fit)["b1"]
      c <- coef(fit)["c1"]

      Topt = -b/(2*a)
      Amax = quadratic.model(Topt, a,b,c)
      d = c-0.5*Amax
      Tlow = (-b - sqrt(b*b-4*a*d))/(2*a)
      Thigh = (-b + sqrt(b*b-4*a*d))/(2*a)
      Tbreadth = abs(Thigh-Tlow)
      CTmin = (-b + sqrt(b*b-4*a*c))/(2*a)
      CTmax = (-b - sqrt(b*b-4*a*c))/(2*a)

      # Compute pseudo-r2
      r_sq = 1-sum(resid(fit)^2)/sum((d_1$A-mean(d_1$A))^2)
    }


    # Bind dataframe
    output.params = bind_rows(output.params,
                              data.frame(curveID = i, Topt, Amax, Tbreadth,
                                         CTmin, CTmax, r_sq))

  }

  return(output.params)
}
#
#
# ###############
# ###############
#
# # Run 'read_6800_and_plot_curves.R' until dataframe 'at.all' is created, or load at.all.csv
# # at.all <- read.csv("../clean_data/at.all.csv") # use this if file is stored locally
# at.all.sub <- at.all %>% select(A, Tleaf, curveID)
#
# # Fit curves based on Arrhenius and Sharpe-Schoolfield models
# fit.info <- fit_curves_michaletz_2021(at.all.sub, x="Tleaf", y="A", T_ref=10)
# fit.info$curveID <- as.character(fit.info$curveID)
# dev.off()
#
# # Fit curves based on quadratic equations
# fits <- fit.quadratic(at.all.sub)
# fits <- subset(fits, Tbreadth < 100)
# fits$curveID <- as.character(fits$curveID)
# head(fits)
#
# # Add meta data back into the fit parameters
# at.meta <- read.csv("../raw_data/meta/A_T_meta.csv")
# at.meta$curveID <- as.character(at.meta$curveID)
# at.meta$site <- as.character(at.meta$site)
#
# # Join Licor curve data with meta-data
# at.join <- left_join(fit.info, at.meta, by="curveID")
# at.join <- left_join(at.join, fits, by="curveID")
#
# # write.csv(at.join, "../clean_data/at.join.csv")
#
#
# # # Join Licor curve data with manual curve quality information
# #
# # at.hit <- read.csv("A_T_curve_hitlist.csv")
# # at.hit <-  at.hit %>% select(curveID, quality)
# # at.hit$curveID <- as.character(at.hit$curveID)
# #
# # at.join2 <- left_join(at.join, at.hit, by="curveID")
# # at.join2 <- subset(at.join, species != "Agrostis")
# #
# # #write.csv(at.join2, "../clean_data/at.join2.csv")
#
#
#
# # # Plots
# # at.join <- read.csv("../clean_data/at.join.csv")
# # at.join$site <- as.character(at.join$site)
# #
# # # subset all the curves for those we deemed good (3)
# # at.join <- subset(at.join, quality == 3)
# #
# #
# # ggplot(data = at.join, aes(x = site, y = Tbreadth, fill = species)) +
# #   geom_dotplot(binaxis="y", stackdir="center", position=position_dodge()) +
# #   theme_bw()
# #
# # ggplot(data = at.join, aes(x = Elevation.masl, y = T_opt, colour = species)) +
# #   geom_point() +
# #   geom_smooth(method="lm", se=FALSE) +
# #   theme_bw()
#

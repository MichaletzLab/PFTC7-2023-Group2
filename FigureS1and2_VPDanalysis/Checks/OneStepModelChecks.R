#This is to confirm the two step method gets the same results as if we run the fits simultaneously.
#This is to confirm the two step method gets the same results as if we run the fits simultaneously.

presence_absence <- dat %>%
  distinct(site, species) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = species, values_from = presence, values_fill = list(presence = 0))

#Three balanced panels possible:
#subset1: Achillea m, Alchemilla a, and vaccinium @ site 6,7,8,9
#subset2: ecklonis and pallidum @ site 1,4
#subset3: pallidum and senecio @ site 2,5

#Make the three balanced panel datasets:
# Recode site and species variables
subset1 <- dat %>%
  filter(country == "Norway" & species != "agrostis capillaris") %>%
  mutate(
    site = factor(site, levels = c(6, 7, 8, 9)),
    species = factor(species, levels = c("achillea millefolium", "alchemilla alpina", "vaccinium vitis-idaea"))
  )


subset2 <- dat %>%
  filter(country == "SAfrica" & 
           site %in% c(1, 4) & 
           species %in% c("helichrysum ecklonis", "helichrysum pallidum"))%>%
  mutate(site = factor(site, levels = c(1,4)))


subset3 <- dat %>%
  filter(country == "SAfrica" & 
           site %in% c(2, 5) & 
           species %in% c("senecio glaberrimus", "helichrysum pallidum"))%>%
  mutate(site = factor(site, levels = c(2,5)))

#Run some test linear regressions to see what I am working with:
summary(lmer(photo~poly(tleaf,2)+tleaf:species+(1|species)+species:elevation+elevation, data=subset1))
summary(lmer(photo~poly(tleaf,2)+tleaf:species+(1|species)+species:elevation+elevation, data=subset2))
summary(lmer(photo~poly(tleaf,2)+tleaf:species+(1|species)+species:elevation+elevation, data=subset3))

#Make dummy columns for species and site in all three subsets

#Subset1:
site_dummies6789 <- model.matrix(~ site, data = subset1)
species_dummies6789 <- model.matrix(~ species, data = subset1)
colnames(site_dummies6789) <- paste0("Elev", seq_along(colnames(site_dummies6789)))
colnames(species_dummies6789) <- paste0("x", seq_along(colnames(species_dummies6789)))
subset1 <- cbind(subset1, site_dummies6789, species_dummies6789)
subset1<-subset1%>%
  mutate(x1 = as.numeric(x1),
         x2 = as.numeric(x2),
         Elev1 = as.numeric(Elev1),
         Elev2 = as.numeric(Elev2),
         Elev3 = as.numeric(Elev3))
#Subset2:
site_dummies14 <- model.matrix(~ site, data = subset2)
species_dummies14 <- model.matrix(~ species, data = subset2)
colnames(site_dummies14) <- paste0("Elev", seq_along(colnames(site_dummies14)))
colnames(species_dummies14) <- paste0("x", seq_along(colnames(species_dummies14)))
subset2 <- cbind(subset2, site_dummies14, species_dummies14)
subset2<-subset2%>%
  mutate(x1 = as.numeric(x1),
         Elev1 = as.numeric(Elev1))

#Subset3:
site_dummies25 <- model.matrix(~ site, data = subset3)
species_dummies25 <- model.matrix(~ species, data = subset3)
colnames(site_dummies25) <- paste0("Elev", seq_along(colnames(site_dummies25)))
colnames(species_dummies25) <- paste0("x", seq_along(colnames(species_dummies25)))
subset3 <- cbind(subset3, site_dummies25, species_dummies25)
subset3<-subset3%>%
  mutate(x1 = as.numeric(x1),
         Elev1 = as.numeric(Elev1))

#Redefine the schoolfield function to have fixed effects for species and elevation
schoolfield_fixed_3s4e <- function(temp, J_ref, topt0, topt1, topt2, toptelev0, toptelev1, toptelev2, toptelev3, E0, E1, E2, Eelev0, Eelev1, Eelev2, Eelev3, Ed0, Ed1, Ed2, Edelev0, Edelev1, Edelev2, Edelev3, x1, x2, Elev1, Elev2, Elev3,T_ref = 25) {
  k <- 8.62e-5
  temp <- temp + 273.15  # Convert to Kelvin
  T_ref <- T_ref + 273.15  # Reference temperature (Kelvin)
  topt0 <- topt0 + 273.15
  topt1<- topt1 + 273.15
  topt2<- topt2 + 273.15
  toptelev0<- toptelev0 + 273.15
  toptelev1<-toptelev1 + 273.15
  toptelev2<- toptelev2 + 273.15
  toptelev3<- toptelev3 + 273.15
  
  # Calculate optimized parameters
  T_opt <- topt0 + topt1 * x1 + topt2 * x2 + toptelev0 + toptelev1 * Elev1 + toptelev2 * Elev2 + toptelev3 * Elev3
  E <- E0 + E1 * x1 + E2 * x2 + Eelev0 + Eelev1 * Elev1 + Eelev2 * Elev2 + Eelev3 * Elev3
  E_D <- Ed0 + Ed1 * x1 + Ed2 * x2 + Edelev0 + Edelev1 * Elev1 + Edelev2 * Elev2 + Edelev3 * Elev3
  
  # Calculate and return Schoolfield model
  return(J_ref * exp(E * (1 / (k * T_ref) - 1 / (k * temp))) / (1 + E / (E_D - E) * exp((E_D / k) * (1 / T_opt - 1 / temp))))
}

schoolfield_fixed_2s2e <- function(temp, J_ref, x1, Elev1, topt0,topt1,toptelev0,toptelev1,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) {
  k = 8.62e-5
  temp <- temp + 273.15  # Convert to Kelvin
  T_ref <- 25 + 273.15  # Reference temperature (Kelvin)
  topt0<- topt0 + 273.15
  topt1<- topt1 + 273.15
  toptelev0<- toptelev0 + 273.15
  toptelev1<- toptelev1 + 273.15
  T_opt <- topt0+topt1*x1 + toptelev0+toptelev1*Elev1
  E <- E0+E1*x1 + Eelev0+Eelev1*Elev1
  E_D <- Ed0+Ed1*x1 + Edelev0+Edelev1*Elev1
  return(J_ref * exp(E * (1/(k*T_ref) - 1/(k*temp))) / 
           (1 + E/(E_D - E) * exp((E_D/k) * (1/T_opt - 1/temp))))
}

#Try modifying the schoolfield nls fit:

fit_mod34_schoolfield <- function(Data, x = "tleaf", x1="x1", x2="x2", Elev1="Elev1", Elev2="Elev2",Elev3="Elev3", y = "photo", T_ref = 25, fitted_models = TRUE) {
  models <- c("schoolfield")
  Data$T_ref <- T_ref
  
  # Standardize names of x and y for use throughout
  if (x != "tleaf") {
    Data$tleaf <- Data[[x]]  # Rename X as tleaf
    Data[[x]] <- NULL        # Get rid of originals
  }
  if (y != "photo") {
    Data$photo <- Data[[y]]  # Rename Y as A
    Data[[y]] <- NULL    # Get rid of originals
  }
  
  Data <- subset(Data, photo > 0)
  
  # Check for missing values
  if (anyNA(Data$tleaf) || anyNA(Data$photo)) {
    stop("Missing values detected in the input data.")
  }
  # Check lengths
  if (length(Data$tleaf) != length(Data$photo)) {
    stop("The lengths of tleaf and photo columns do not match.")
  }
  
  # Initialize a list to store the results
  schoolfield.results.data <- list()
  
  ######################
  # Curve fitting      #
  ######################
  
  print("Fitting Schoolfield model to the entire dataset")
  
  plot <- ggplot(Data, aes(x = tleaf, y = photo)) +
    geom_point() +
    theme_classic()
  
  # Fit Schoolfield function using nls_multstart
  fit <- tryCatch({
    nls.multstart::nls_multstart(photo ~ schoolfield_fixed_3s4e(temp = tleaf, J_ref, topt0, topt1, topt2, toptelev0, toptelev1, toptelev2, toptelev3, E0, E1, E2, Eelev0, Eelev1, Eelev2, Eelev3, Ed0, Ed1, Ed2, Edelev0, Edelev1, Edelev2, Edelev3, x1, x2, Elev1, Elev2, Elev3),
                                 data = Data,
                                 iter = 200,
                                 start_lower = c(J_ref = 0, topt0 = 10, topt1 = 0, topt2 = 0, toptelev0 = 10, toptelev1 = 0, toptelev2 = 0, toptelev3 = 0, E0 = 0, E1 = 0, E2 = 0, Eelev0 = 0, Eelev1 = 0, Eelev2 = 0, Eelev3 = 0, Ed0 = 0, Ed1 = 0, Ed2 = 0, Edelev0 = 0, Edelev1 = 0, Edelev2 = 0, Edelev3 = 0),
                                 start_upper = c(J_ref = 20, topt0 = 40, topt1 = 10, topt2 = 10, toptelev0 = 40, toptelev1 = 10, toptelev2 = 10, toptelev3 = 10, E0 = 10, E1 = 2, E2 = 2, Eelev0 = 2, Eelev1 = 2, Eelev2 = 2, Eelev3 = 2, Ed0 = 15, Ed1 = 5, Ed2 = 5, Edelev0 = 15, Edelev1 = 5, Edelev2 = 5, Edelev3 = 5),
                                 supp_errors = 'Y',
                                 na.action = na.omit,
                                 lower = c(J_ref = 0, topt0 = 10, topt1 = 0, topt2 = 0, toptelev0 = 10, toptelev1 = 0, toptelev2 = 0, toptelev3 = 0, E0 = 0, E1 = 0, E2 = 0, Eelev0 = 0, Eelev1 = 0, Eelev2 = 0, Eelev3 = 0, Ed0 = 0, Ed1 = 0, Ed2 = 0, Edelev0 = 0, Edelev1 = 0, Edelev2 = 0, Edelev3 = 0))
  }, error = function(e) {
    print("Fitting failed for the entire dataset")
    print(e)
    return(NULL)
  })
  
  if (is.null(fit)) {
    print("Fit is NULL, returning empty dataframe")
    return(data.frame())
  }
  
  print("Fit was successful")
  
  # Store predictions for schoolfield
  predicted_schoolfield <- fitted(fit)
  
  # Calculate R-squared for schoolfield
  r_sq <- 1 - sum(resid(fit)^2) / sum((Data$photo - mean(Data$photo))^2)
  
  # Extract parameters estimates and SE
  params <- coef(fit)
  params_SE <- summary(fit)$coefficients[,'Std. Error']
  AIC <- AIC(fit)
  
  # Predict A values for the entire range of tleaf
  predicted_photo <- predict(fit, newdata = data.frame(tleaf = Data$tleaf))
  
  
  schoolfield.results.data <- data.frame(
    J_ref = params["J_ref"], J_ref_SE = params_SE["J_ref"],
    topt0 = params["topt0"], topt0_SE = params_SE["topt0"],
    topt1 = params["topt1"], topt1_SE = params_SE["topt1"],
    topt2 = params["topt2"], topt2_SE = params_SE["topt2"],
    toptelev0 = params["toptelev0"], toptelev0_SE = params_SE["toptelev0"],
    toptelev1 = params["toptelev1"], toptelev1_SE = params_SE["toptelev1"],
    toptelev2 = params["toptelev2"], toptelev2_SE = params_SE["toptelev2"],
    toptelev3 = params["toptelev3"], toptelev3_SE = params_SE["toptelev3"],
    E0 = params["E0"], E0_SE = params_SE["E0"],
    E1 = params["E1"], E1_SE = params_SE["E1"],
    E2 = params["E2"], E2_SE = params_SE["E2"],
    Eelev0 = params["Eelev0"], Eelev0_SE = params_SE["Eelev0"],
    Eelev1 = params["Eelev1"], Eelev1_SE = params_SE["Eelev1"],
    Eelev2 = params["Eelev2"], Eelev2_SE = params_SE["Eelev2"],
    Eelev3 = params["Eelev3"], Eelev3_SE = params_SE["Eelev3"],
    Ed0 = params["Ed0"], Ed0_SE = params_SE["Ed0"],
    Ed1 = params["Ed1"], Ed1_SE = params_SE["Ed1"],
    Ed2 = params["Ed2"], Ed2_SE = params_SE["Ed2"],
    Edelev0 = params["Edelev0"], Edelev0_SE = params_SE["Edelev0"],
    Edelev1 = params["Edelev1"], Edelev1_SE = params_SE["Edelev1"],
    Edelev2 = params["Edelev2"], Edelev2_SE = params_SE["Edelev2"],
    Edelev3 = params["Edelev3"], Edelev3_SE = params_SE["Edelev3"],
    AIC = AIC, r_sq = r_sq, breadth = breadth
  )
  
  print(schoolfield.results.data)  # Print the result dataframe for debugging
  
  # Return the results data frame
  return(schoolfield.results.data)
}

SchoolfieldMODIFIED34 <- fit_mod34_schoolfield(subset1, x = "tleaf", x1="x1", x2="x2", Elev1="Elev1", Elev2="Elev2",Elev3="Elev3", y = "photo", T_ref = 25)

























#Create the two new basis functions:
schoolfield_basis3x4 <- function(data, x_var, J_ref, x1, x2, Elev1, Elev2, Elev3, topt0, topt1, topt2, toptelev0, toptelev1, toptelev2, toptelev3, E0, E1, E2, Eelev0, Eelev1, Eelev2, Eelev3, Ed0, Ed1, Ed2, Edelev0, Edelev1, Edelev2, Edelev3) {
  x <- data[[x_var]]
  basis <- schoolfield_fixed_3s4e(x, J_ref, x1, x2, Elev1, Elev2, Elev3, topt0, topt1, topt2, toptelev0, toptelev1, toptelev2, toptelev3, E0, E1, E2, Eelev0, Eelev1, Eelev2, Eelev3, Ed0, Ed1, Ed2, Edelev0, Edelev1, Edelev2, Edelev3) 
  return(basis)
}
schoolfield_basis2x2 <- function(data, x_var,J_ref, x1, Elev1, topt0,topt1,toptelev0,toptelev1,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) {
  x <- data[[x_var]]
  basis <- schoolfield_fixed_2s2e(x, J_ref, x1, Elev1, topt0,topt1,toptelev0,toptelev1,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) 
  return(basis)
}



## run subset1: ##
fit_schoolfield_gam_fixed34 <- function(data, x_var, y_var, start_params) {
  # Objective function to minimize (negative log-likelihood)
  obj_fun <- function(params) {
    J_ref <- params[1]
    topt0 <- params[2]
    topt1 <- params[3]
    topt2 <- params[4]
    toptelev0 <- params[5]
    toptelev1 <- params[6]
    toptelev2 <- params[7]
    toptelev3 <- params[8]
    E0 <- params[9]
    E1 <- params[10]
    E2 <- params[11]
    Eelev0 <- params[12]
    Eelev1 <- params[13]
    Eelev2 <- params[14]
    Eelev3 <- params[15]
    Ed0 <- params[16]
    Ed1 <- params[17]
    Ed2 <- params[18]
    Edelev0 <- params[19]
    Edelev1 <- params[20]
    Edelev2 <- params[21]
    Edelev3 <- params[22]
    
    # Calculate the basis for the GAM
    data$schoolfield_base <- schoolfield_fixed_3s4e(data[[x_var]], J_ref, topt0, topt1, topt2, toptelev0, toptelev1, toptelev2, toptelev3, E0, E1, E2, Eelev0, Eelev1, Eelev2, Eelev3, Ed0, Ed1, Ed2, Edelev0, Edelev1, Edelev2, Edelev3, data$x1, data$x2, data$Elev1, data$Elev2, data$Elev3)
    
    # Fit the GAM model
    gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
    
    # Print the negative log-likelihood
    nll <- -logLik(gam_fit)
    print(paste("Negative Log-Likelihood:", nll))
    print(params)
    
    return(nll)
  }
  
  # Optimize parameters with different methods for better convergence
  opt_res <- optim(start_params, obj_fun, method = "L-BFGS-B", lower = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), upper = c(40, 40, 40, 40, 40, 40, 40, 40, 5, 5, 5, 5, 5, 5, 5, 15, 15, 15, 15, 15, 15, 15))
  
  # Print optimization results for debugging
  print(opt_res)
  
  # Check if the optimization converged
  if (opt_res$convergence != 0) {
    stop("Optimization did not converge")
  }
  
  # Extract optimized parameters
  opt_params <- opt_res$par
  
  # Fit GAM model with optimized parameters
  data$schoolfield_base <- schoolfield_fixed_3s4e(data[[x_var]], opt_params[1], opt_params[2], opt_params[3], opt_params[4], opt_params[5], opt_params[6], opt_params[7], opt_params[8], opt_params[9], opt_params[10], opt_params[11], opt_params[12], opt_params[13], opt_params[14], opt_params[15], opt_params[16], opt_params[17], opt_params[18], opt_params[19], opt_params[20], opt_params[21], opt_params[22], data$x1, data$x2, data$Elev1, data$Elev2, data$Elev3)
  
  gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
  
  return(list(gam_fit = gam_fit, opt_params = opt_params))
}

# Example usage
start_params <- c(J_ref = 5, topt0 = 20, topt1 = 1, topt2 = 1, toptelev0 = 22, toptelev1 = 1, toptelev2 = 1, toptelev3 = 1, E0 = 0.5, E1 = 0.25, E2 = 0.25, Eelev0 = 0.5, Eelev1 = 0.25, Eelev2 = 0.25, Eelev3 = 0.25, Ed0 = 1.5, Ed1 = 0.5, Ed2 = 0.5, Edelev0 = 1.5, Edelev1 = 0.5, Edelev2 = 0.5, Edelev3 = 0.5) 
fit_fixedSchoolfield34 <- fit_schoolfield_gam_fixed34(subset1, "tleaf", "photo", start_params = start_params)
summary(fit_fixedSchoolfield34)

opt_params <- fit_fixedSchoolfield34$opt_params
topt0 <- opt_params["topt0"]
topt1 <- opt_params["topt1"]
topt2 <- opt_params["topt2"]











## run subset2 and subset3: ##
fit_schoolfield_gam_fixed22 <- function(data, x_var, y_var, start_params) {
  #data <- na.omit(data[, c(x_var, y_var, "x1", "Elev1")])
  # Objective function to minimize (negative log-likelihood)
  obj_fun <- function(params) {
    J_ref <- params[1]
    topt0 <- params[2]
    topt1 <- params[3]
    toptelev0 <- params[4]
    toptelev1 <- params[5]
    E0 <- params[6]
    E1 <- params[7]
    Eelev0 <- params[8]
    Eelev1 <- params[9]
    Ed0 <- params[10]
    Ed1 <- params[11]
    Edelev0 <- params[12]
    Edelev1 <- params[13]
    
    # Calculate the basis for the GAM
    data$schoolfield_base <- schoolfield_fixed_2s2e(data[[x_var]], J_ref, topt0, topt1, toptelev0, toptelev1, E0, E1, Eelev0, Eelev1, Ed0, Ed1, Edelev0, Edelev1, data$x1, data$Elev1)
    
    # Fit the GAM model
    gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
    
    # Print the negative log-likelihood
    nll <- -logLik(gam_fit)
    print(paste("Negative Log-Likelihood:", nll))
    print(params)
    
    return(nll)
  }
  
  # Optimize parameters with different methods for better convergence
  opt_res <- optim(start_params, obj_fun, method = "L-BFGS-B", lower = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), upper = c(40, 40, 40, 40, 40, 5, 5, 5, 5, 15, 15, 15, 15))
  
  # Print optimization results for debugging
  print(opt_res)
  
  # Check if the optimization converged
  if (opt_res$convergence != 0) {
    stop("Optimization did not converge")
  }
  
  # Extract optimized parameters
  opt_params <- opt_res$par
  
  # Fit GAM model with optimized parameters
  data$schoolfield_base <- schoolfield_fixed_2s2e(data[[x_var]], opt_params[1], opt_params[2], opt_params[3], opt_params[4], opt_params[5], opt_params[6], opt_params[7], opt_params[8], opt_params[9], opt_params[10], opt_params[11], opt_params[12], opt_params[13], data$x1, data$Elev1)
  
  gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
  
  return(list(gam_fit = gam_fit, opt_params = opt_params))
}

# Example usage
start_params <- c(J_ref = 5, topt0 = 20, topt1 = 1, toptelev0 = 20, toptelev1 = 1, E0 = 1, E1 = 0.01, Eelev0 = 0.5, Eelev1 = 0.25, Ed0 = 1.5, Ed1 = 0.5, Edelev0 = 1.5, Edelev1 = 0.1) 

fit_fixedSchoolfield14 <- fit_schoolfield_gam_fixed22(subset2, "tleaf", "photo", start_params = start_params)
fit_fixedSchoolfield25 <- fit_schoolfield_gam_fixed22(subset3, "tleaf", "photo", start_params = start_params)

summary(fit_fixedSchoolfield22)

opt_params <- fit_fixedSchoolfield14$opt_params
(topt0 <- opt_params["topt0"])
(topt1 <- opt_params["topt1"])
(toptelev0 <- opt_params["toptelev0"])
(toptelev1 <- opt_params["toptelev1"])

(E0 <- opt_params["E0"])
(E1 <- opt_params["E1"])
(Eelev0 <- opt_params["Eelev0"])
(Eelev1 <- opt_params["Eelev1"])

(Ed0 <- opt_params["Ed0"])
(Ed1 <- opt_params["Ed1"])
(Edelev0 <- opt_params["Edelev0"])
(Edelev1 <- opt_params["Edelev1"])









################### This is working sort of for topt with no intercept term. However, topt2 in both categories is just the starting value... ##############
schoolfield_fixed_2s2e <- function(temp, J_ref, x1,x2, Elev1, Elev2,topt1,topt2,toptelev1,toptelev2,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) {
  k = 8.62e-5
  temp <- temp + 273.15  # Convert to Kelvin
  T_ref <- 25 + 273.15  # Reference temperature (Kelvin)
  topt2<- topt2 + 273.15
  topt1<- topt1 + 273.15
  toptelev2<- toptelev2 + 273.15
  toptelev1<- toptelev1 + 273.15
  T_opt <- topt1*x1 + topt2*x2 + toptelev1*Elev1 + toptelev2*Elev2
  E <- E0+E1*x1 + Eelev0+Eelev1*Elev1
  E_D <- Ed0+Ed1*x1 + Edelev0+Edelev1*Elev1
  return(J_ref * exp(E * (1/(k*T_ref) - 1/(k*temp))) / 
           (1 + E/(E_D - E) * exp((E_D/k) * (1/T_opt - 1/temp))))
}
schoolfield_basis2x2 <- function(data, x_var,J_ref, x1,x2, Elev1, Elev2,topt1,topt2,toptelev1,toptelev2,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) {
  x <- data[[x_var]]
  basis <- schoolfield_fixed_2s2e(x, J_ref, x1,x2, Elev1, Elev2,topt0,topt1,toptelev0,toptelev1,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) 
  return(basis)
}
fit_schoolfield_gam_fixed22 <- function(data, x_var, y_var, start_params) {
  #data <- na.omit(data[, c(x_var, y_var, "x1", "Elev1")])
  # Objective function to minimize (negative log-likelihood)
  obj_fun <- function(params) {
    J_ref <- params[1]
    topt1 <- params[2]
    topt2 <- params[3]
    toptelev1 <- params[4]
    toptelev2 <- params[5]
    E0 <- params[6]
    E1 <- params[7]
    Eelev0 <- params[8]
    Eelev1 <- params[9]
    Ed0 <- params[10]
    Ed1 <- params[11]
    Edelev0 <- params[12]
    Edelev1 <- params[13]
    
    # Calculate the basis for the GAM
    data$schoolfield_base <- schoolfield_fixed_2s2e(data[[x_var]], J_ref, topt1, topt2, toptelev1, toptelev2, E0, E1, Eelev0, Eelev1, Ed0, Ed1, Edelev0, Edelev1, data$x1, data$x2, data$Elev1,data$Elev2)
    
    # Fit the GAM model
    gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base-1")), data = data, method = "REML")
    
    # Print the negative log-likelihood
    nll <- -logLik(gam_fit)
    print(paste("Negative Log-Likelihood:", nll))
    print(params)
    
    return(nll)
  }
  
  # Optimize parameters with different methods for better convergence
  opt_res <- optim(start_params, obj_fun, method = "L-BFGS-B", lower = c(0, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0), upper = c(10, 30, 30, 30, 30, 5, 5, 5, 5, 15, 15, 15, 15))
  
  # Print optimization results for debugging
  print(opt_res)
  
  # Check if the optimization converged
  if (opt_res$convergence != 0) {
    stop("Optimization did not converge")
  }
  
  # Extract optimized parameters
  opt_params <- opt_res$par
  
  # Fit GAM model with optimized parameters
  data$schoolfield_base <- schoolfield_fixed_2s2e(data[[x_var]], opt_params[1], opt_params[2], opt_params[3], opt_params[4], opt_params[5], opt_params[6], opt_params[7], opt_params[8], opt_params[9], opt_params[10], opt_params[11], opt_params[12], opt_params[13], data$x1,  data$x2, data$Elev1,data$Elev2)
  
  gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base-1")), data = data, method = "REML")
  
  return(list(gam_fit = gam_fit, opt_params = opt_params))
}

# Example usage
start_params <- c(J_ref = 3, topt1 = 20, topt2 = 15, toptelev1 = 20, toptelev2 = 15, E0 = 0.1, E1 = 0.01, Eelev0 = 0.5, Eelev1 = 1, Ed0 = 1.5, Ed1 = 0.5, Edelev0 = 1.5, Edelev1 = 0.1) 

fit_fixedSchoolfield14 <- fit_schoolfield_gam_fixed22(subset2, "tleaf", "photo", start_params = start_params)
fit_fixedSchoolfield25 <- fit_schoolfield_gam_fixed22(subset3, "tleaf", "photo", start_params = start_params)

summary(fit_fixedSchoolfield22)

opt_params <- fit_fixedSchoolfield25$opt_params
(J_ref <- opt_params["J_ref"])
(topt1 <- opt_params["topt1"])
(topt2 <- opt_params["topt2"])
(toptelev1 <- opt_params["toptelev1"])
(toptelev2 <- opt_params["toptelev2"])

(E0 <- opt_params["E0"])
(E1 <- opt_params["E1"])
(Eelev0 <- opt_params["Eelev0"])
(Eelev1 <- opt_params["Eelev1"])

(Ed0 <- opt_params["Ed0"])
(Ed1 <- opt_params["Ed1"])
(Edelev0 <- opt_params["Edelev0"])
(Edelev1 <- opt_params["Edelev1"])








################### Not working... ##############
schoolfield_fixed_2s2e <- function(temp, J_ref, x1,x2, Elev1, Elev2,topt1,topt2,toptelev1,toptelev2,E,E_D, T_ref=25) {
  k = 8.62e-5
  temp <- temp + 273.15  # Convert to Kelvin
  T_ref <- T_ref + 273.15  # Reference temperature (Kelvin)
  topt1<- topt1 + 273.15
  topt2<- topt2 + 273.15
  toptelev1<- toptelev1 + 273.15
  toptelev2<- toptelev2 + 273.15
  T_opt <- topt1*x1 + topt2*x2 + toptelev1*Elev1 + toptelev2*Elev2
  #E <- E0+E1*x1 + Eelev0+Eelev1*Elev1
  #E_D <- Ed0+Ed1*x1 + Edelev0+Edelev1*Elev1
  return(J_ref * exp(E * (1/(k*T_ref) - 1/(k*temp))) / 
           (1 + E/(E_D - E) * exp((E_D/k) * (1/T_opt - 1/temp))))
}
schoolfield_basis2x2 <- function(data, x_var,J_ref, x1,x2, Elev1, Elev2,topt1,topt2,toptelev1,toptelev2,E,E_D) {
  x <- data[[x_var]]
  basis <- schoolfield_fixed_2s2e(x, J_ref, x1,x2, Elev1, Elev2,topt1,topt2,toptelev1,toptelev2,E,E_D) 
  return(basis)
}
fit_schoolfield_gam_fixed22 <- function(data, x_var, y_var, start_params) {
  #data <- na.omit(data[, c(x_var, y_var, "x1", "Elev1")])
  # Objective function to minimize (negative log-likelihood)
  obj_fun <- function(params) {
    J_ref <- params[1]
    topt1 <- params[2]
    topt2 <- params[3]
    toptelev1 <- params[4]
    toptelev2 <- params[5]
    E <- params[6]
    E_D <- params[7]
    
    # Calculate the basis for the GAM
    data$schoolfield_base <- schoolfield_fixed_2s2e(data[[x_var]], J_ref, topt1, topt2, toptelev1, toptelev2, E, E_D, data$x1, data$x2, data$Elev1,data$Elev2)
    
    # Fit the GAM model
    gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
    
    # Print the negative log-likelihood
    nll <- -logLik(gam_fit)
    print(paste("Negative Log-Likelihood:", nll))
    print(params)
    
    return(nll)
  }
  
  # Optimize parameters with different methods for better convergence
  opt_res <- optim(start_params, obj_fun, method = "L-BFGS-B", lower = c(0, 15, 15, 15, 15, 0.01, 0.01), upper = c(10, 30, 30, 30, 30, 5, 10))
  
  # Print optimization results for debugging
  print(opt_res)
  
  # Check if the optimization converged
  if (opt_res$convergence != 0) {
    stop("Optimization did not converge")
  }
  
  # Extract optimized parameters
  opt_params <- opt_res$par
  
  # Fit GAM model with optimized parameters
  data$schoolfield_base <- schoolfield_fixed_2s2e(data[[x_var]], opt_params[1], opt_params[2], opt_params[3], opt_params[4], opt_params[5], opt_params[6], opt_params[7], data$x1,  data$x2, data$Elev1,data$Elev2)
  
  gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
  
  return(list(gam_fit = gam_fit, opt_params = opt_params))
}

# Example usage
start_params <- c(J_ref = 3, topt1 = 20, topt2 = 15, toptelev1 = 20, toptelev2 = 15, E = 0.5, E_D = 1.5) 

fit_fixedSchoolfield14 <- fit_schoolfield_gam_fixed22(subset2, "tleaf", "photo", start_params = start_params)
fit_fixedSchoolfield25 <- fit_schoolfield_gam_fixed22(subset3, "tleaf", "photo", start_params = start_params)

summary(fit_fixedSchoolfield22)

opt_params <- fit_fixedSchoolfield14$opt_params
(J_ref <- opt_params["J_ref"])
(topt1 <- opt_params["topt1"])
(topt2 <- opt_params["topt2"])
(toptelev1 <- opt_params["toptelev1"])
(toptelev2 <- opt_params["toptelev2"])

(E0 <- opt_params["E0"])
(E1 <- opt_params["E1"])
(Eelev0 <- opt_params["Eelev0"])
(Eelev1 <- opt_params["Eelev1"])

(Ed0 <- opt_params["Ed0"])
(Ed1 <- opt_params["Ed1"])
(Edelev0 <- opt_params["Edelev0"])
(Edelev1 <- opt_params["Edelev1"])









######This is a newer version where pretty much nothing works::
presence_absence <- dat %>%
  distinct(site, species) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = species, values_from = presence, values_fill = list(presence = 0))

#Three balanced panels possible:
#subset1: Achillea m, Alchemilla a, and vaccinium @ site 6,7,8,9
#subset2: ecklonis and pallidum @ site 1,4
#subset3: pallidum and senecio @ site 2,5

#Make the three balanced panel datasets:
# Recode site and species variables
subset1 <- dat %>%
  filter(country == "Norway" & species != "agrostis capillaris") %>%
  mutate(
    site = factor(site, levels = c(6, 7, 8, 9)),
    species = factor(species, levels = c("achillea millefolium", "alchemilla alpina", "vaccinium vitis-idaea"))
  )


subset2 <- dat %>%
  filter(country == "SAfrica" & 
           site %in% c(1, 4) & 
           species %in% c("helichrysum ecklonis", "helichrysum pallidum"))%>%
  mutate(site = factor(site, levels = c(1,4)))


subset3 <- dat %>%
  filter(country == "SAfrica" & 
           site %in% c(2, 5) & 
           species %in% c("senecio glaberrimus", "helichrysum pallidum"))%>%
  mutate(site = factor(site, levels = c(2,5)))

#Run some test linear regressions to see what I am working with:
summary(lmer(photo~poly(tleaf,2)+tleaf:species+(1|species)+species:elevation+elevation, data=subset1))
summary(lmer(photo~poly(tleaf,2)+tleaf:species+(1|species)+species:elevation+elevation, data=subset2))
summary(lmer(photo~poly(tleaf,2)+tleaf:species+(1|species)+species:elevation+elevation, data=subset3))

#Make dummy columns for species and site in all three subsets

#Subset1:
site_dummies6789 <- model.matrix(~ site, data = subset1)
species_dummies6789 <- model.matrix(~ species, data = subset1)
colnames(site_dummies6789) <- paste0("Elev", seq_along(colnames(site_dummies6789)))
colnames(species_dummies6789) <- paste0("x", seq_along(colnames(species_dummies6789)))
subset1 <- cbind(subset1, site_dummies6789, species_dummies6789)
subset1<-subset1%>%
  mutate(x1 = as.numeric(x1),
         x2 = as.numeric(x2),
         Elev1 = as.numeric(Elev1),
         Elev2 = as.numeric(Elev2),
         Elev3 = as.numeric(Elev3))
#Subset2:
site_dummies14 <- model.matrix(~ site, data = subset2)
species_dummies14 <- model.matrix(~ species, data = subset2)
colnames(site_dummies14) <- paste0("Elev", seq_along(colnames(site_dummies14)))
colnames(species_dummies14) <- paste0("x", seq_along(colnames(species_dummies14)))
subset2 <- cbind(subset2, site_dummies14, species_dummies14)
subset2<-subset2%>%
  mutate(x1 = as.numeric(x1),
         Elev1 = as.numeric(Elev1))

#Subset3:
site_dummies25 <- model.matrix(~ site, data = subset3)
species_dummies25 <- model.matrix(~ species, data = subset3)
colnames(site_dummies25) <- paste0("Elev", seq_along(colnames(site_dummies25)))
colnames(species_dummies25) <- paste0("x", seq_along(colnames(species_dummies25)))
subset3 <- cbind(subset3, site_dummies25, species_dummies25)
subset3<-subset3%>%
  mutate(x1 = as.numeric(x1),
         Elev1 = as.numeric(Elev1))

#Redefine the schoolfield function to have fixed effects for species and elevation
schoolfield_fixed_3s4e <- function(temp, J_ref, topt0, topt1, topt2, toptelev0, toptelev1, toptelev2, toptelev3, E0, E1, E2, Eelev0, Eelev1, Eelev2, Eelev3, Ed0, Ed1, Ed2, Edelev0, Edelev1, Edelev2, Edelev3, x1, x2, Elev1, Elev2, Elev3,T_ref = 25) {
  k <- 8.62e-5
  temp <- temp + 273.15  # Convert to Kelvin
  T_ref <- T_ref + 273.15  # Reference temperature (Kelvin)
  topt0 <- topt0 + 273.15
  topt1<- topt1 + 273.15
  topt2<- topt2 + 273.15
  toptelev0<- toptelev0 + 273.15
  toptelev1<-toptelev1 + 273.15
  toptelev2<- toptelev2 + 273.15
  toptelev3<- toptelev3 + 273.15
  
  # Calculate optimized parameters
  T_opt <- topt0 + topt1 * x1 + topt2 * x2 + toptelev0 + toptelev1 * Elev1 + toptelev2 * Elev2 + toptelev3 * Elev3
  E <- E0 + E1 * x1 + E2 * x2 + Eelev0 + Eelev1 * Elev1 + Eelev2 * Elev2 + Eelev3 * Elev3
  E_D <- Ed0 + Ed1 * x1 + Ed2 * x2 + Edelev0 + Edelev1 * Elev1 + Edelev2 * Elev2 + Edelev3 * Elev3
  
  # Calculate and return Schoolfield model
  return(J_ref * exp(E * (1 / (k * T_ref) - 1 / (k * temp))) / (1 + E / (E_D - E) * exp((E_D / k) * (1 / T_opt - 1 / temp))))
}

schoolfield_fixed_2s2e <- function(temp, J_ref, x1,x2, Elev1, Elev2,topt1,topt2,toptelev1,toptelev2,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) {
  k = 8.62e-5
  temp <- temp + 273.15  # Convert to Kelvin
  T_ref <- 25 + 273.15  # Reference temperature (Kelvin)
  topt2<- topt2 + 273.15
  topt1<- topt1 + 273.15
  toptelev2<- toptelev2 + 273.15
  toptelev1<- toptelev1 + 273.15
  T_opt <- topt1*x1 + topt2*x2 + toptelev1*Elev1 + toptelev2*Elev2
  E <- E0+E1*x1 + Eelev0+Eelev1*Elev1
  E_D <- Ed0+Ed1*x1 + Edelev0+Edelev1*Elev1
  return(J_ref * exp(E * (1/(k*T_ref) - 1/(k*temp))) / 
           (1 + E/(E_D - E) * exp((E_D/k) * (1/T_opt - 1/temp))))
}

schoolfield_basis2x2 <- function(data, x_var,J_ref, x1,x2, Elev1, Elev2,topt1,topt2,toptelev1,toptelev2,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) {
  x <- data[[x_var]]
  basis <- schoolfield_fixed_2s2e(x, J_ref, x1,x2, Elev1, Elev2,topt0,topt1,toptelev0,toptelev1,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) 
  return(basis)
}

fit_schoolfield_gam_fixed22 <- function(data, x_var, x1,x2,Elev1,Elev2,y_var, start_params) {
  
  # Objective function to minimize (negative log-likelihood)
  obj_fun <- function(params) {
    J_ref <- params[1]
    topt1 <- params[2]
    topt2 <- params[3]
    toptelev1 <- params[4]
    toptelev2 <- params[5]
    E0 <- params[6]
    E1 <- params[7]
    Eelev0 <- params[8]
    Eelev1 <- params[9]
    Ed0 <- params[10]
    Ed1 <- params[11]
    Edelev0 <- params[12]
    Edelev1 <- params[13]
    
    # Calculate the basis for the GAM
    data$schoolfield_base <- schoolfield_fixed_2s2e(data[[x_var]], J_ref, topt1, topt2, toptelev1, toptelev2, E0, E1, Eelev0, Eelev1, Ed0, Ed1, Edelev0, Edelev1, data$x1,data$x2, data$Elev1,data$Elev2)
    
    # Fit the GAM model
    gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
    
    # Print the negative log-likelihood
    nll <- -logLik(gam_fit)
    print(paste("Negative Log-Likelihood:", nll))
    print(params)
    
    return(nll)
  }
  
  # Optimize parameters with different methods for better convergence
  opt_res <- optim(start_params, obj_fun, method = "L-BFGS-B", lower = c(0, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0), upper = c(10, 30, 30, 30, 30, 5, 5, 5, 5, 10, 10, 10, 10), control = list(fnscale = -1))
  
  # Print optimization results for debugging
  print(opt_res)
  
  # Check if the optimization converged
  if (opt_res$convergence != 0) {
    stop("Optimization did not converge")
  }
  
  # Extract optimized parameters
  opt_params <- opt_res$par
  
  # Fit GAM model with optimized parameters
  data$schoolfield_base <- schoolfield_fixed_2s2e(data[[x_var]], opt_params[1], opt_params[2], opt_params[3], opt_params[4], opt_params[5], opt_params[6], opt_params[7], opt_params[8], opt_params[9], opt_params[10], opt_params[11], opt_params[12], opt_params[13], data$x1, data$x2, data$Elev1,data$Elev2)
  
  gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
  
  return(list(gam_fit = gam_fit, opt_params = opt_params))
}

# Example usage
start_params <- c(J_ref = 5, topt1 = 20, topt2 = 20, toptelev1 = 20, toptelev2 = 20, E0 = 0, E1 = 0, Eelev0 = 0.5, Eelev1 = 0.25, Ed0 = 1.5, Ed1 = 0.5, Edelev0 = 0, Edelev1 = 0) 
fit_fixedSchoolfield22 <- fit_schoolfield_gam_fixed22(subset1, "tleaf", "x1","x2","Elev1","Elev2","photo", start_params = start_params)

























#Try modifying the schoolfield nls fit:

fit_mod34_schoolfield <- function(Data, x = "tleaf", x1="x1", x2="x2", Elev1="Elev1", Elev2="Elev2",Elev3="Elev3", y = "photo", T_ref = 25, fitted_models = TRUE) {
  models <- c("schoolfield")
  Data$T_ref <- T_ref
  
  # Standardize names of x and y for use throughout
  if (x != "tleaf") {
    Data$tleaf <- Data[[x]]  # Rename X as tleaf
    Data[[x]] <- NULL        # Get rid of originals
  }
  if (y != "photo") {
    Data$photo <- Data[[y]]  # Rename Y as A
    Data[[y]] <- NULL    # Get rid of originals
  }
  
  Data <- subset(Data, photo > 0)
  
  # Check for missing values
  if (anyNA(Data$tleaf) || anyNA(Data$photo)) {
    stop("Missing values detected in the input data.")
  }
  # Check lengths
  if (length(Data$tleaf) != length(Data$photo)) {
    stop("The lengths of tleaf and photo columns do not match.")
  }
  
  # Initialize a list to store the results
  schoolfield.results.data <- list()
  
  ######################
  # Curve fitting      #
  ######################
  
  print("Fitting Schoolfield model to the entire dataset")
  
  plot <- ggplot(Data, aes(x = tleaf, y = photo)) +
    geom_point() +
    theme_classic()
  
  # Fit Schoolfield function using nls_multstart
  fit <- tryCatch({
    nls.multstart::nls_multstart(photo ~ schoolfield_fixed_3s4e(temp = tleaf, J_ref, topt0, topt1, topt2, toptelev0, toptelev1, toptelev2, toptelev3, E0, E1, E2, Eelev0, Eelev1, Eelev2, Eelev3, Ed0, Ed1, Ed2, Edelev0, Edelev1, Edelev2, Edelev3, x1, x2, Elev1, Elev2, Elev3),
                                 data = Data,
                                 iter = 200,
                                 start_lower = c(J_ref = 0, topt0 = 10, topt1 = 0, topt2 = 0, toptelev0 = 10, toptelev1 = 0, toptelev2 = 0, toptelev3 = 0, E0 = 0, E1 = 0, E2 = 0, Eelev0 = 0, Eelev1 = 0, Eelev2 = 0, Eelev3 = 0, Ed0 = 0, Ed1 = 0, Ed2 = 0, Edelev0 = 0, Edelev1 = 0, Edelev2 = 0, Edelev3 = 0),
                                 start_upper = c(J_ref = 20, topt0 = 40, topt1 = 10, topt2 = 10, toptelev0 = 40, toptelev1 = 10, toptelev2 = 10, toptelev3 = 10, E0 = 10, E1 = 2, E2 = 2, Eelev0 = 2, Eelev1 = 2, Eelev2 = 2, Eelev3 = 2, Ed0 = 15, Ed1 = 5, Ed2 = 5, Edelev0 = 15, Edelev1 = 5, Edelev2 = 5, Edelev3 = 5),
                                 supp_errors = 'Y',
                                 na.action = na.omit,
                                 lower = c(J_ref = 0, topt0 = 10, topt1 = 0, topt2 = 0, toptelev0 = 10, toptelev1 = 0, toptelev2 = 0, toptelev3 = 0, E0 = 0, E1 = 0, E2 = 0, Eelev0 = 0, Eelev1 = 0, Eelev2 = 0, Eelev3 = 0, Ed0 = 0, Ed1 = 0, Ed2 = 0, Edelev0 = 0, Edelev1 = 0, Edelev2 = 0, Edelev3 = 0))
  }, error = function(e) {
    print("Fitting failed for the entire dataset")
    print(e)
    return(NULL)
  })
  
  if (is.null(fit)) {
    print("Fit is NULL, returning empty dataframe")
    return(data.frame())
  }
  
  print("Fit was successful")
  
  # Store predictions for schoolfield
  predicted_schoolfield <- fitted(fit)
  
  # Calculate R-squared for schoolfield
  r_sq <- 1 - sum(resid(fit)^2) / sum((Data$photo - mean(Data$photo))^2)
  
  # Extract parameters estimates and SE
  params <- coef(fit)
  params_SE <- summary(fit)$coefficients[,'Std. Error']
  AIC <- AIC(fit)
  
  # Predict A values for the entire range of tleaf
  predicted_photo <- predict(fit, newdata = data.frame(tleaf = Data$tleaf))
  
  
  schoolfield.results.data <- data.frame(
    J_ref = params["J_ref"], J_ref_SE = params_SE["J_ref"],
    topt0 = params["topt0"], topt0_SE = params_SE["topt0"],
    topt1 = params["topt1"], topt1_SE = params_SE["topt1"],
    topt2 = params["topt2"], topt2_SE = params_SE["topt2"],
    toptelev0 = params["toptelev0"], toptelev0_SE = params_SE["toptelev0"],
    toptelev1 = params["toptelev1"], toptelev1_SE = params_SE["toptelev1"],
    toptelev2 = params["toptelev2"], toptelev2_SE = params_SE["toptelev2"],
    toptelev3 = params["toptelev3"], toptelev3_SE = params_SE["toptelev3"],
    E0 = params["E0"], E0_SE = params_SE["E0"],
    E1 = params["E1"], E1_SE = params_SE["E1"],
    E2 = params["E2"], E2_SE = params_SE["E2"],
    Eelev0 = params["Eelev0"], Eelev0_SE = params_SE["Eelev0"],
    Eelev1 = params["Eelev1"], Eelev1_SE = params_SE["Eelev1"],
    Eelev2 = params["Eelev2"], Eelev2_SE = params_SE["Eelev2"],
    Eelev3 = params["Eelev3"], Eelev3_SE = params_SE["Eelev3"],
    Ed0 = params["Ed0"], Ed0_SE = params_SE["Ed0"],
    Ed1 = params["Ed1"], Ed1_SE = params_SE["Ed1"],
    Ed2 = params["Ed2"], Ed2_SE = params_SE["Ed2"],
    Edelev0 = params["Edelev0"], Edelev0_SE = params_SE["Edelev0"],
    Edelev1 = params["Edelev1"], Edelev1_SE = params_SE["Edelev1"],
    Edelev2 = params["Edelev2"], Edelev2_SE = params_SE["Edelev2"],
    Edelev3 = params["Edelev3"], Edelev3_SE = params_SE["Edelev3"],
    AIC = AIC, r_sq = r_sq, breadth = breadth
  )
  
  print(schoolfield.results.data)  # Print the result dataframe for debugging
  
  # Return the results data frame
  return(schoolfield.results.data)
}

SchoolfieldMODIFIED34 <- fit_mod34_schoolfield(subset1, x = "tleaf", x1="x1", x2="x2", Elev1="Elev1", Elev2="Elev2",Elev3="Elev3", y = "photo", T_ref = 25)























#########TRYING THIS::::::##########

#Create the two new basis functions:
schoolfield_basis3x4 <- function(data, x_var, J_ref, x1, x2, Elev1, Elev2, Elev3, topt0, topt1, topt2, toptelev0, toptelev1, toptelev2, toptelev3, E0, E1, E2, Eelev0, Eelev1, Eelev2, Eelev3, Ed0, Ed1, Ed2, Edelev0, Edelev1, Edelev2, Edelev3) {
  x <- data[[x_var]]
  basis <- schoolfield_fixed_3s4e(x, J_ref, x1, x2, Elev1, Elev2, Elev3, topt0, topt1, topt2, toptelev0, toptelev1, toptelev2, toptelev3, E0, E1, E2, Eelev0, Eelev1, Eelev2, Eelev3, Ed0, Ed1, Ed2, Edelev0, Edelev1, Edelev2, Edelev3) 
  return(basis)
}
schoolfield_basis2x2 <- function(data, x_var,J_ref, x1, Elev1, topt0,topt1,toptelev0,toptelev1,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) {
  x <- data[[x_var]]
  basis <- schoolfield_fixed_2s2e(x, J_ref, x1, Elev1, topt0,topt1,toptelev0,toptelev1,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) 
  return(basis)
}



## run subset1: ##
fit_schoolfield_gam_fixed34 <- function(data, x_var, y_var, start_params) {
  # Objective function to minimize (negative log-likelihood)
  obj_fun <- function(params) {
    J_ref <- params[1]
    topt0 <- params[2]
    topt1 <- params[3]
    topt2 <- params[4]
    toptelev0 <- params[5]
    toptelev1 <- params[6]
    toptelev2 <- params[7]
    toptelev3 <- params[8]
    E0 <- params[9]
    E1 <- params[10]
    E2 <- params[11]
    Eelev0 <- params[12]
    Eelev1 <- params[13]
    Eelev2 <- params[14]
    Eelev3 <- params[15]
    Ed0 <- params[16]
    Ed1 <- params[17]
    Ed2 <- params[18]
    Edelev0 <- params[19]
    Edelev1 <- params[20]
    Edelev2 <- params[21]
    Edelev3 <- params[22]
    
    # Calculate the basis for the GAM
    data$schoolfield_base <- schoolfield_fixed_3s4e(data[[x_var]], J_ref, topt0, topt1, topt2, toptelev0, toptelev1, toptelev2, toptelev3, E0, E1, E2, Eelev0, Eelev1, Eelev2, Eelev3, Ed0, Ed1, Ed2, Edelev0, Edelev1, Edelev2, Edelev3, data$x1, data$x2, data$Elev1, data$Elev2, data$Elev3)
    
    # Fit the GAM model
    gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
    
    # Print the negative log-likelihood
    nll <- -logLik(gam_fit)
    print(paste("Negative Log-Likelihood:", nll))
    print(params)
    
    return(nll)
  }
  
  # Optimize parameters with different methods for better convergence
  opt_res <- optim(start_params, obj_fun, method = "L-BFGS-B", lower = c(0, 15, 15, 15, 15, 15, 15, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), upper = c(10, 30, 30, 30, 30, 30, 30, 30, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5))
  
  # Print optimization results for debugging
  print(opt_res)
  
  # Check if the optimization converged
  if (opt_res$convergence != 0) {
    stop("Optimization did not converge")
  }
  
  # Extract optimized parameters
  opt_params <- opt_res$par
  
  # Fit GAM model with optimized parameters
  data$schoolfield_base <- schoolfield_fixed_3s4e(data[[x_var]], opt_params[1], opt_params[2], opt_params[3], opt_params[4], opt_params[5], opt_params[6], opt_params[7], opt_params[8], opt_params[9], opt_params[10], opt_params[11], opt_params[12], opt_params[13], opt_params[14], opt_params[15], opt_params[16], opt_params[17], opt_params[18], opt_params[19], opt_params[20], opt_params[21], opt_params[22], data$x1, data$x2, data$Elev1, data$Elev2, data$Elev3)
  
  gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
  
  return(list(gam_fit = gam_fit, opt_params = opt_params))
}

# Example usage
start_params <- c(J_ref = 5, topt0 = 25, topt1 = 1, topt2 = 1, toptelev0 = 25, toptelev1 = 1, toptelev2 = 1, toptelev3 = 1, E0 = 1, E1 = 0.01, E2 = 001, Eelev0 = 0.5, Eelev1 = 0.25, Eelev2 = 0.25, Eelev3 = 0.25, Ed0 = 1.5, Ed1 = 0.5, Ed2 = 0.5, Edelev0 = 0.5, Edelev1 = 0.05, Edelev2 = 0.05, Edelev3 = 0.05) 
fit_fixedSchoolfield34 <- fit_schoolfield_gam_fixed34(subset1, "tleaf", "photo", start_params = start_params)
summary(fit_fixedSchoolfield34)

opt_params <- fit_fixedSchoolfield34$opt_params
(topt0 <- opt_params["topt0"])
(topt1 <- opt_params["topt1"])
(topt2 <- opt_params["topt2"])
(toptelev0 <- opt_params["toptelev0"])
(toptelev1 <- opt_params["toptelev1"])
(toptelev2 <- opt_params["toptelev2"])
(toptelev3 <- opt_params["toptelev3"])

























schoolfield_fixed_2s2e <- function(temp, J_ref, x1, x2, Elev1, Elev2, topt1, topt2, toptelev1, toptelev2, E, E_D, T_ref = 25) {
  k = 8.62e-5
  temp <- temp + 273.15  # Convert to Kelvin
  T_ref <- T_ref + 273.15  # Reference temperature (Kelvin)
  topt2 <- topt2 + 273.15
  topt1 <- topt1 + 273.15
  toptelev2 <- toptelev2 + 273.15
  toptelev1 <- toptelev1 + 273.15
  T_opt <- topt1 * x1 + topt2 * x2 + toptelev1 * Elev1 + toptelev2 * Elev2
  return(J_ref * exp(E * (1 / (k * T_ref) - 1 / (k * temp))) / 
           (1 + E / (E_D - E) * exp((E_D / k) * (1 / T_opt - 1 / temp))))
}

schoolfield_basis2x2 <- function(data, x_var, J_ref, x1, x2, Elev1, Elev2, topt1, topt2, toptelev1, toptelev2, E, E_D) {
  x <- data[[x_var]]
  basis <- schoolfield_fixed_2s2e(x, J_ref, x1, x2, Elev1, Elev2, topt1, topt2, toptelev1, toptelev2, E, E_D)
  return(basis)
}





## run subset2 and subset3: ##
fit_schoolfield_gam_fixed22 <- function(data, x_var, y_var, start_params) {
  # Check if the required columns are present in the data
  required_columns <- c(x_var, y_var, "x1", "x2", "Elev1", "Elev2")
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(paste("Missing columns in the data:", paste(missing_columns, collapse = ", ")))
  }
  
  # Remove rows with NA values in the relevant columns
  data <- na.omit(data[, ..required_columns])
  
  # Objective function to minimize (negative log-likelihood)
  obj_fun <- function(params) {
    J_ref <- params[1]
    topt1 <- params[2]
    topt2 <- params[3]
    toptelev1 <- params[4]
    toptelev2 <- params[5]
    E <- params[6]
    E_D <- params[7]
    
    # Calculate the basis for the GAM
    data$schoolfield_base <- schoolfield_fixed_2s2e(data[[x_var]], J_ref, data$x1, data$x2, data$Elev1, data$Elev2, topt1, topt2, toptelev1, toptelev2, E, E_D)
    
    # Fit the GAM model
    gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
    
    # Calculate the negative log-likelihood
    nll <- -logLik(gam_fit)
    
    return(nll)
  }
  
  # Define lower and upper bounds
  lower_bounds <- c(1, 15, 15, 15, 15, 0, 0)
  upper_bounds <- c(8, 30, 30, 30, 30, 5, 5)
  
  # Run the optimization using DEoptim
  opt_res <- DEoptim(obj_fun, lower = lower_bounds, upper = upper_bounds, DEoptim.control(NP = 70, itermax = 100))
  
  # Extract optimized parameters
  opt_params <- opt_res$optim$bestmem
  
  # Print optimization results for debugging
  print(opt_res)
  
  # Check if the optimization converged
  if (!opt_res$optim$convergence) {
    stop("Optimization did not converge")
  }
  
  # Fit GAM model with optimized parameters
  data$schoolfield_base <- schoolfield_fixed_2s2e(data[[x_var]], opt_params[1], data$x1, data$x2, data$Elev1, data$Elev2, opt_params[2], opt_params[3], opt_params[4], opt_params[5], opt_params[6], opt_params[7])
  
  gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
  
  return(list(gam_fit = gam_fit, opt_params = opt_params))
}

# Example usage
start_params <- c(J_ref = 5, topt1 = 20, topt2 = 20, toptelev1 = 20, toptelev2 = 20, E = 0.5, E_D = 1.5) 

fit_fixedSchoolfield14 <- fit_schoolfield_gam_fixed22(subset2, "tleaf", "photo", start_params = start_params)



fit_fixedSchoolfield25 <- fit_schoolfield_gam_fixed22(subset3, "tleaf", "photo", start_params = start_params)

summary(fit_fixedSchoolfield22)

opt_params <- fit_fixedSchoolfield14$opt_params
(J_ref <- opt_params["J_ref"])
(topt1 <- opt_params["topt1"])
(topt2 <- opt_params["topt2"])
(toptelev2 <- opt_params["toptelev2"])
(toptelev1 <- opt_params["toptelev1"])

(E <- opt_params["E"])
(E_D <- opt_params["E_D"])









#This is getting the error: Error in if (sum(edf2) > sum(edf1)) { : missing value where TRUE/FALSE needed

schoolfield_fixed_2s2e <- function(temp, J_ref, x1,x2, Elev1, Elev2,topt0,topt1,toptelev0,toptelev1,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) {
  k = 8.62e-5
  temp <- temp + 273.15  # Convert to Kelvin
  T_ref <- 25 + 273.15  # Reference temperature (Kelvin)
  topt2<- topt2 + 273.15
  topt1<- topt1 + 273.15
  toptelev2<- toptelev2 + 273.15
  toptelev1<- toptelev1 + 273.15
  T_opt <- topt1*x1 + topt2*x2 + toptelev1*Elev1 + toptelev2*Elev2
  E <- E0+E1*x1 + Eelev0+Eelev1*Elev1
  E_D <- Ed0+Ed1*x1 + Edelev0+Edelev1*Elev1
  return(J_ref * exp(E * (1/(k*T_ref) - 1/(k*temp))) / 
           (1 + E/(E_D - E) * exp((E_D/k) * (1/T_opt - 1/temp))))
}

schoolfield_basis2x2 <- function(data, x_var,J_ref, x1,x2, Elev1, Elev2,topt1,topt2,toptelev1,toptelev2,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) {
  x <- data[[x_var]]
  basis <- schoolfield_fixed_2s2e(x, J_ref, x1,x2, Elev1, Elev2,topt0,topt1,toptelev0,toptelev1,E0,E1,Eelev0,Eelev1,Ed0,Ed1,Edelev0,Edelev1) 
  return(basis)
}

fit_schoolfield_gam_fixed22 <- function(data, x_var, y_var, start_params) {
  
  # Objective function to minimize (negative log-likelihood)
  obj_fun <- function(params) {
    J_ref <- params[1]
    topt1 <- params[2]
    topt2 <- params[3]
    toptelev1 <- params[4]
    toptelev2 <- params[5]
    E0 <- params[6]
    E1 <- params[7]
    Eelev0 <- params[8]
    Eelev1 <- params[9]
    Ed0 <- params[10]
    Ed1 <- params[11]
    Edelev0 <- params[12]
    Edelev1 <- params[13]
    
    # Calculate the basis for the GAM
    data$schoolfield_base <- schoolfield_fixed_2s2e(data[[x_var]], J_ref, topt1, topt2, toptelev1, toptelev2, E0, E1, Eelev0, Eelev1, Ed0, Ed1, Edelev0, Edelev1, data$x1,data$x2, data$Elev1,data$Elev2)
    
    # Fit the GAM model
    gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
    
    # Print the negative log-likelihood
    nll <- -logLik(gam_fit)
    print(paste("Negative Log-Likelihood:", nll))
    print(params)
    
    return(nll)
  }
  
  # Optimize parameters with different methods for better convergence
  opt_res <- optim(start_params, obj_fun, method = "L-BFGS-B", lower = c(0, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0), upper = c(10, 30, 30, 30, 30, 5, 5, 5, 5, 10, 10, 10, 10), control = list(fnscale = -1))
  
  # Print optimization results for debugging
  print(opt_res)
  
  # Check if the optimization converged
  if (opt_res$convergence != 0) {
    stop("Optimization did not converge")
  }
  
  # Extract optimized parameters
  opt_params <- opt_res$par
  
  # Fit GAM model with optimized parameters
  data$schoolfield_base <- schoolfield_fixed_2s2e(data[[x_var]], opt_params[1], opt_params[2], opt_params[3], opt_params[4], opt_params[5], opt_params[6], opt_params[7], opt_params[8], opt_params[9], opt_params[10], opt_params[11], opt_params[12], opt_params[13], data$x1, data$x2, data$Elev1,data$Elev1)
  
  gam_fit <- gam(as.formula(paste(y_var, "~ schoolfield_base")), data = data, method = "REML")
  
  return(list(gam_fit = gam_fit, opt_params = opt_params))
}

# Example usage
start_params <- c(J_ref = 5, topt0 = 20, topt1 = 1, toptelev0 = 10, toptelev1 = 0, E0 = 0, E1 = 0, Eelev0 = 0.5, Eelev1 = 0.25, Ed0 = 1.5, Ed1 = 0.5, Edelev0 = 0, Edelev1 = 0) 
fit_fixedSchoolfield22 <- fit_schoolfield_gam_fixed22(subset1, "tleaf", "photo", start_params = start_params)


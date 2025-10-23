# Purpose: Run a gam on all data - adding structure as needed
# Plots: Photosynthesis vs. tleaf showing species intercepts
        #Photosynthesis vs. tleaf faceted by species and colored for elevation 
# Dependencies: configure.datfile.R

#######################################
########### Model Fitting ###########
GAM_SS0 = gam(photo~s(tleaf,k=5)+s(elevation, k=5), data=dats,method = "REML")
summary(GAM_SS0)

GAM_SS2 = gam(photo~s(tleaf,k=5)+s(elevation, k=5) +factor(curveid), data=dats,method = "REML")
summary(GAM_SS2)

datsss = dat%>%
  filter(country=="SAfrica")

GAM_SS3 = gam(photo~s(tleaf,k=5)+s(elevation, k=5), data=datsss,method = "REML")
summary(GAM_SS3)

# Run a mixed effects model to account for a random nested effect of curveid and species
LMER_SS = lmer(photo ~ poly(tleaf,2) + elevation+ elevation:species+ (1|species), data = dats)
summary(LMER_SS)

LMER_SS=lmer(elevation ~ (1|species), data=dats)
r_squared <- r.squaredGLMM(LMER_SS)
r_squared
#######################################
########### Prediction Data ###########
dats$species <- factor(dats$species)
species_combinations <- data.frame(species = levels(dats$species))

# Generate prediction data with all combinations of tleaf, elevation, and species
pred_data <- expand.grid(
  tleaf = seq(min(dats$tleaf), max(dats$tleaf), length.out = 100),
  elevation = seq(min(dats$elevation), max(dats$elevation), length.out = 100),
  species = species_combinations$species
)

# Ensure species in pred_data matches the levels of species in dats
pred_data$species <- factor(pred_data$species, levels = levels(dats$species))

# Predict using the GAM model
pred_data$predicted <- predict(GAM_SS3, newdata = pred_data, type = "response")



#######################################
###########     Plotting   ###########
# Plot with species as separate lines on one plot - shows intercept difference between species
(p <- ggplot(pred_data, aes(x = tleaf, y = predicted, color = species)) +
  geom_smooth(size = 1, aes(fill = species), alpha = 0.2) +
  labs(
    x = "Leaf temperature (˚C)",
    y = "Predicted photosynthesis") +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  ))


#Now plot faceting by species and coloring for elevation
(p <- ggplot(pred_data, aes(x = tleaf, y = predicted, color = elevation)) +
  geom_point(size = 2, alpha = 0.6) +  # Scatter plot with color gradient for elevation
  labs(
    x = "Leaf Temperature (˚C)",
    y = "Predicted Photosynthesis"
  ) +
  scale_color_gradient(low = "red", high = "blue") +  # Color gradient for elevation
  theme_minimal()+
  facet_wrap(~species, scales = "free_y"))
# Purpose: Run a gam on all data - adding structure as needed
# Plots: Photosynthesis vs. tleaf showing species intercepts
        #Photosynthesis vs. tleaf faceted by species and colored for elevation 

########### Model Fitting ###########
GAM_SS0 = gam(photo~s(tleaf,k=5)+s(elevation, k=5), data=dats,method = "REML")
summary(GAM_SS0)

#GAM_SS1 = gam(photo~predicted+s(elevation, k=5), data=dats,method = "REML") # I don't think this is a good idea since predicted becomes a linear term and gets rid of important error
#summary(GAM_SS1)

GAM_SS2 = gam(photo~s(tleaf,k=5)+s(elevation, k=5) +factor(curveid), data=dats,method = "REML") # I don't think we should be adding a leaf-level effect
summary(GAM_SS2)

GAM_SS3 = gam(photo~s(tleaf,k=5)+s(elevation, k=5) + factor(species), data=dats,method = "REML")
summary(GAM_SS3)

########### Prediction Data ###########
# Convert species to a factor with appropriate levels
dats$species <- factor(dats$species)

# Check the levels again
levels(dats$species)
# Create a data frame with all unique combinations of species
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

# Plotting
p <- ggplot(pred_data, aes(x = tleaf, y = predicted, color = species)) +
  geom_smooth(size = 1, aes(fill = species), alpha = 0.2) +  # Set fill for error bands
  labs(
    x = "Leaf temperature (˚C)",
    y = "Predicted photosynthesis") +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

print(p)



p <- ggplot(pred_data, aes(x = tleaf, y = predicted, color = elevation)) +
  geom_point(size = 2, alpha = 0.6) +  # Scatter plot with color gradient for elevation
  labs(
    x = "Leaf Temperature (˚C)",
    y = "Predicted Photosynthesis"
  ) +
  scale_color_gradient(low = "red", high = "blue") +  # Color gradient for elevation
  theme_minimal()+
  facet_wrap(~species, scales = "free_y")

print(p)


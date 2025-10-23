# Fit Schoolfield GAM all in one method
ssmod.tleaf <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[dat$mod_group == sel_group, ]
  fit_results <- fit_schoolfield_gam("tleaf", "photo", T_ref = 25, start_params = start_params, data = tmp_dat)
  return(list(mod_group = sel_group, gam_fit = fit_results$gam_fit))
})


# Extract smooth estimates for Schoolfield GAM
ssmod_sm_tleaf <- lapply(ssmod.tleaf, function(x) {
  smooth_estimates(x$gam_fit, unconditional = TRUE,smooth = 's(tleaf)', partial_match = TRUE) %>%
    mutate(mod = "sst", R2 = summary(x$gam_fit)$r.sq)
}) %>% rbindlist()



####### Liu's Method #############
# Define the Schoolfield model function
schoolfield <- function(temp, J_ref, E, E_D, T_opt, T_ref = 25) {
  k <- 8.62e-5            # Boltzmann's constant (eV K-1)
  temp <- temp + 273.15    # Convert to Kelvin
  T_ref <- T_ref + 273.15  # Convert to Kelvin
  T_opt <- T_opt + 273.15  # Convert to Kelvin
  return(J_ref * exp(E * (1/(k*T_ref) - 1/(k*temp))) / (1 + E/(E_D - E) * exp((E_D/k) * (1/T_opt - 1/temp))))
}

# Merge dat with schoolfield fits based on curveid
dats <- dat %>%
  left_join(schoolfield.fit, by = "curveid") %>%
  mutate(predicted = schoolfield(tleaf, J_ref, E, E_D, T_opt.s))

vec_mod_group <- unique(dats$mod_group)

# Define the function to fit the GAM model and extract smooth estimates for Liu's method
fit_gam_with_schoolfield_liu <- function(data, x_var, y_var, predicted) {
  # Ensure data is clean and complete
  if (anyNA(names(data)) || anyNA(data[[x_var]]) || anyNA(data[[y_var]])) {
    stop("Data contains NA values in critical columns.")
  }
  
  # Fit the GAM model with s.tleaf and predicted
  formula <- as.formula(paste(y_var, "~ s(", x_var, ") + s(predicted)"))
  gam_fit <- gam(formula, data = data, method = "REML")
  
  # Extract smooth estimates
  smooth_estimates <- predict(gam_fit, type = "terms", se.fit = TRUE)
  smooth_df <- data.frame(.estimate = smooth_estimates$fit[, "s(tleaf)"], .se = smooth_estimates$se.fit[, "s(tleaf)"])
  smooth_df$tleaf <- data[[x_var]]
  
  return(list(gam_fit = gam_fit, smooth_df = smooth_df))
}

# Iterate over mod groups, fit the model, and extract smooth estimates
ssmod.tleaf_liu <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dats[dats$mod_group == sel_group, ]
  fit_results <- tryCatch(
    fit_gam_with_schoolfield_liu(tmp_dat, "tleaf", "photo","predicted"),
    error = function(e) NULL
  )
  if (!is.null(fit_results)) {
    smooth_df <- fit_results$smooth_df %>%
      mutate(mod = "liu", R2 = summary(fit_results$gam_fit)$r.sq)
    return(smooth_df)
  } else {
    return(NULL)
  }
}) %>% bind_rows()

# Combine all data
ssmod_data <- rbindlist(list(ssmod_sm_tleaf, ssmod.tleaf_liu), fill = TRUE)

# Ensure data is ordered by tleaf
ssmod_data <- ssmod_data %>%
  group_by(mod) %>%
  arrange(tleaf)

## Plot without faceting
p_liu.comp <- ssmod_data %>%
  ggplot(aes(x = tleaf, y = .estimate, color = mod, fill = mod)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "ts"), se = TRUE, size = 1, alpha = 0.25) +
  geom_hline(aes(yintercept = 0), col = 'grey40', size = 0.5, linetype = 1) +
  scale_color_manual(values = c("sst" = "orange", "liu" = "purple4"),
                     labels = c("sst" = "GAM Schoolfield + Tleaf", "liu" = "nlsSchoolfield + Tleaf")) +
  scale_fill_manual(values = c("sst" = "orange", "liu" = "purple4"),
                    labels = c("sst" = "GAM Schoolfield + Tleaf", "liu" = "nlsSchoolfield + Tleaf")) +
  scale_x_continuous(breaks = seq(5, 40, 5), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-10, 15))+
  labs(color = 'GAM', fill = 'GAM', 
       y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m^-2, s^-1, ")")), 
       x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(legend.position = 'bottom') #+
  #geom_vline(data = avg_weib_topt, aes(xintercept = avg_weib_topt), color = 'grey', linetype = 2, show.legend = FALSE, size = 1)

# Save the Plot
ggsave(p_liu.comp, 
       filename = paste0("pftc7_vpd_analysis/figures/Liu_comparison_plot_", Sys.Date(), ".png"),
       device = "png",
       width = 40,
       height = 40,
       units = 'cm',
       scale = 0.8,
       dpi = 600)

#### g1 plot ####

# --- Function to fit g1 per curveID ---
fit_g1 <- function(df){
  # filter out bad/missing data
  uso_data <- df %>% filter(!is.na(gsw), !is.na(A), !is.na(Ca), !is.na(VPDleaf), gsw > 0, A > 0, Ca > 0, VPDleaf > 0)
  
  # require at least 10 points
  if(nrow(uso_data) < 10) return(NA)
  
  # fit the USO model
  fit <- tryCatch(
    nls(gsw ~ g0 + 1.6 * (1 + g1 / sqrt(VPDleaf)) * (A / Ca),
        data = uso_data,
        start = list(g0 = 0.01, g1 = 4)),
    error = function(e) NA
  )
  
  # return g1 coefficient
  if(is.na(fit)[1]) return(NA)
  coef(fit)["g1"]
}

# --- Calculate g1 for each curveID ---
g1_curve <- raw.dat %>%
  group_by(curveID) %>%
  summarise(
    g1 = fit_g1(cur_data()),
    Elevation = first(Elevation),   # grab Elevation for this curve
    .groups = "drop"
  )

# --- Join mean air temperature by elevation ---
g1_curve <- g1_curve %>%
  left_join(
    air_by_elev_source %>% select(Elevation, mean_air),
    by = "Elevation"
  )

# --- Check result ---
head(g1_curve)

# --- Plot g1 vs mean air temperature ---
ggplot(g1_curve, aes(x = mean_air, y = g1)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 5),
    se = TRUE,         # show confidence interval
    color = "blue",
    fill = "lightblue"
  ) +
  labs(
    x = "Mean Air Temperature (°C)",
    y = expression(g[1])
  ) +
  theme_classic()
lm_g1 <- lm(g1 ~ mean_air, data = g1_curve)
summary(lm_g1) # not different than 0





#### g1 vs Soil Moisture ####

# --- Calculate g1 per curveID (same as before) ---
g1_curve_moist <- raw.dat %>%
  group_by(curveID) %>%
  summarise(
    g1 = fit_g1(cur_data()),
    Elevation = first(Elevation),   # keep Elevation for join
    .groups = "drop"
  )

# --- Join mean soil moisture by elevation ---
g1_curve_moist <- g1_curve_moist %>%
  left_join(
    moist_by_elev_source %>% 
      group_by(Elevation) %>%
      summarise(mean_moist_pct = mean(mean_moist_pct, na.rm = TRUE), .groups = "drop"),
    by = "Elevation"
  )

# --- Check result ---
head(g1_curve_moist)

# --- Plot g1 vs mean soil moisture ---
ggplot(g1_curve_moist, aes(x = mean_moist_pct, y = g1)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 5),
    se = TRUE,
    color = "forestgreen",
    fill = "lightgreen"
  ) +
  labs(
    x = "Mean Soil Moisture (%)",
    y = expression(g[1])
  ) +
  theme_classic()

# --- Optional: linear model for reference ---
lm_g1_moist <- lm(g1 ~ mean_moist_pct, data = g1_curve_moist)
summary(lm_g1_moist)

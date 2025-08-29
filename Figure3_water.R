# Bin soil moisture: -----------------------------------------------------------
gam_A   <- gam(A   ~ moist_bin + s(Tleaf, by = moist_bin, k = 10), data = dat_environ)
gam_E   <- gam(E   ~ moist_bin + s(Tleaf, by = moist_bin, k = 10), data = dat_environ)
gam_gsw <- gam(gsw ~ moist_bin + s(Tleaf, by = moist_bin, k = 10), data = dat_environ)

# --- 1) Ensure moist_bin is factor ---
dat_environ <- dat_environ %>%
  filter(!is.na(moist_bin)) %>%
  mutate(moist_bin = factor(moist_bin, levels = c("Low","High")))

# --- 2) Fit GAMs ---
gam_A   <- gam(A   ~ moist_bin + s(Tleaf, by = moist_bin, k = 10), data = dat_environ)
gam_E   <- gam(E   ~ moist_bin + s(Tleaf, by = moist_bin, k = 10), data = dat_environ)
gam_gsw <- gam(gsw ~ moist_bin + s(Tleaf, by = moist_bin, k = 10), data = dat_environ)

# --- 3) Create prediction grid with standard errors ---
pred_grid <- dat_environ %>%
  group_by(moist_bin) %>%
  summarise(Tleaf = seq(min(Tleaf, na.rm=TRUE),
                        max(Tleaf, na.rm=TRUE),
                        length.out = 100), .groups = "drop")

# Predict with SE for A
pred_A <- predict(gam_A, newdata = pred_grid, se.fit = TRUE)
pred_grid <- pred_grid %>%
  mutate(A_fit = pred_A$fit,
         A_upper = pred_A$fit + pred_A$se.fit,
         A_lower = pred_A$fit - pred_A$se.fit)

# Predict with SE for E
pred_E <- predict(gam_E, newdata = pred_grid, se.fit = TRUE)
pred_grid <- pred_grid %>%
  mutate(E_fit = pred_E$fit,
         E_upper = pred_E$fit + pred_E$se.fit,
         E_lower = pred_E$fit - pred_E$se.fit)

# Predict with SE for gsw
pred_gsw <- predict(gam_gsw, newdata = pred_grid, se.fit = TRUE)
pred_grid <- pred_grid %>%
  mutate(gsw_fit = pred_gsw$fit,
         gsw_upper = pred_gsw$fit + pred_gsw$se.fit,
         gsw_lower = pred_gsw$fit - pred_gsw$se.fit)

# --- 4) Define colors ---
moist_colors <- c("Low" = "lightblue3", "High" = "#3333CC")

# --- 5) Function to plot with ribbon ---
plot_gam <- function(yvar, fitvar, lowervar, uppervar, ylab){
  ggplot() +
    geom_point(data = dat_environ, aes(x = Tleaf, y = .data[[yvar]], color = moist_bin),
               alpha = 0.01) +
    geom_ribbon(data = pred_grid,
                aes(x = Tleaf, ymin = .data[[lowervar]], ymax = .data[[uppervar]], fill = moist_bin),
                alpha = 0.2, inherit.aes = FALSE) +
    geom_line(data = pred_grid, aes(x = Tleaf, y = .data[[fitvar]], color = moist_bin), size = 1.2) +
    scale_color_manual(values = moist_colors) +
    scale_fill_manual(values = moist_colors) +
    labs(x = "Leaf Temperature (°C)", y = ylab, color = "Soil moisture", fill = "Soil moisture") +
    theme_classic()
}

# --- 6) Make individual plots ---
pA    <- plot_gam("A", "A_fit", "A_lower", "A_upper", expression(A~(mu*mol~CO[2]~m^-2~s^-1)))
pE    <- plot_gam("E", "E_fit", "E_lower", "E_upper", expression(E~(mmol~m^-2~s^-1)))
pgsw  <- plot_gam("gsw", "gsw_fit", "gsw_lower", "gsw_upper", expression(g[s][w]~(mol~m^-2~s^-1)))


# --- 3) Function to fit Medlyn g1 ---
fit_g1 <- function(df){
  uso_data <- df %>%
    filter(!is.na(gsw), !is.na(A), !is.na(Ca), !is.na(VPDleaf),
           gsw > 0, A > 0, Ca > 0, VPDleaf > 0)
  
  if(nrow(uso_data) < 10) return(NA)  # safety check
  
  fit <- tryCatch({
    nls(gsw ~ g0 + 1.6 * (1 + g1 / sqrt(VPDleaf)) * (A / Ca),
        data = uso_data,
        start = list(g0 = 0.01, g1 = 4))
  }, error = function(e) NA)
  
  if(is.na(fit)[1]) return(NA)
  coef(fit)["g1"]
}

# --- 4) Calculate g1 per soil moisture bin ---
g1_values <- dat_environ %>%
  group_by(moist_bin) %>%
  summarise(g1 = fit_g1(cur_data()), .groups = "drop")

print(g1_values)

# --- 5) JOIN g1_values to dat_environ and THEN calculate USO_x ---
dat_environ <- dat_environ %>%
  left_join(g1_values, by = "moist_bin")
dat_environ <- dat_environ %>%
  mutate(USO_x = 1.6 * (1 + g1 / sqrt(VPDleaf)) * A / Ca)  # now g1 exists

# --- 6) Build annotation text ---
annot_text <- paste0("g1 Low = ", round(g1_values$g1[g1_values$moist_bin=="Low"], 2),
                     "\n g1 High = ", round(g1_values$g1[g1_values$moist_bin=="High"], 2))

# --- 7) Plot USO model with annotation ---
pUSO <- ggplot(dat_environ, aes(x = USO_x, y = gsw, color = moist_bin)) +
  geom_point(alpha = 0.01) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = moist_colors, name = "Soil moisture") +
  labs(x = expression(1.6~(1+g[1]/sqrt(VPDleaf))*A/C[a]),
       y = expression(g[s][w]~(mol~m^{-2}~s^{-1}))) +
  theme_classic() +
  annotate("text", x = min(dat_environ$USO_x, na.rm = TRUE),
           y = max(dat_environ$gsw, na.rm = TRUE),
           label = annot_text,
           hjust = 0, vjust = 1, size = 4, color = "black")

ggarrange(pA, pE, pgsw, pUSO, nrow = 2, ncol = 2, common.legend = TRUE, labels = c("A","B","C","D"))
#Test if the slopes are the same or not:
lm_USO <- lm(gsw ~ USO_x * moist_bin, data = dat_environ)
summary(lm_USO) #Confirms slopes are the same while intercepts are different.

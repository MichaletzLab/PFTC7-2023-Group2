# --- 1) Fit GAM on full data (with interaction so each air_cat has its own smooth) ---
gam_iWUE_air <- mgcv::gam(I(A/gsw) ~ s(Tleaf, by = air_cat, bs = "cs") + air_cat, 
                          data = dat.airTemp)

# --- 2) Prediction grid across full Tleaf range, by air_cat ---
pred_grid_air <- expand.grid(
  Tleaf = seq(min(dat.airTemp$Tleaf, na.rm = TRUE),
              max(dat.airTemp$Tleaf, na.rm = TRUE),
              length.out = 500),
  air_cat = levels(dat.airTemp$air_cat)
)

pred_air <- predict(gam_iWUE_air, newdata = pred_grid_air, se.fit = TRUE)

pred_grid_air <- pred_grid_air %>%
  mutate(
    y_pred = pred_air$fit,
    se = pred_air$se.fit,
    ymin = y_pred - 1.96 * se,
    ymax = y_pred + 1.96 * se
  )

# --- Find Topt per air_cat ---
topt_by_air <- pred_grid_air %>%
  group_by(air_cat) %>%
  slice_max(order_by = y_pred, n = 1, with_ties = FALSE) %>%
  select(air_cat, Tleaf, y_pred)

# --- Fit overall GAM (no grouping, global curve) ---
gam_iWUE_air_all <- mgcv::gam(I(A/gsw) ~ s(Tleaf, bs = "cs"), data = dat.airTemp)
pred_air_all <- data.frame(Tleaf = seq(min(dat.airTemp$Tleaf, na.rm = TRUE),
                                       max(dat.airTemp$Tleaf, na.rm = TRUE),
                                       length.out = 500))
pred_air_all$y_pred <- predict(gam_iWUE_air_all, newdata = pred_air_all)

overall_topt_air <- pred_air_all[which.max(pred_air_all$y_pred), "Tleaf"]

# --- Define colors for air temperature categories ---
air_colors <- c("Low" = "grey30",
                "Medium" = "lightblue3",
                "High" = "#3333CC")

# --- 3) Filter only for display (y < 400) ---
dat_plot_air <- dat.airTemp %>% filter((A/gsw) < 400)
pred_grid_air_plot <- pred_grid_air %>% filter(y_pred < 400)

# --- 4) Plot with confidence ribbon ---
iWUE.Tleaf.plot.air <- ggplot() +
  geom_point(data = dat_plot_air,
             aes(x = Tleaf, y = A/gsw, color = air_cat),
             alpha = 0.01) +
  geom_ribbon(data = pred_grid_air_plot,
              aes(x = Tleaf, ymin = ymin, ymax = ymax, fill = air_cat),
              alpha = 0.6) +
  geom_line(data = pred_grid_air_plot,
            aes(x = Tleaf, y = y_pred, color = air_cat),
            size = 1.2) +
  geom_vline(xintercept = overall_topt_air, linetype = "dashed", size = 1) +
  geom_vline(data = topt_by_air,
             aes(xintercept = Tleaf, color = air_cat),
             linetype = "dotted", size = 1) +
  labs(
    x = "Leaf temperature (°C)",
    y = expression(iWUE~(mu*mol~CO[2]~mol^{-1}~H[2]*O)),
    color = "Air temperature",
    fill = "Air temperature"
  ) +
  scale_color_manual(values = air_colors) +
  scale_fill_manual(values = air_colors) +
  theme_classic() +
  coord_cartesian(ylim = c(0, 400))

iWUE.Tleaf.plot.air

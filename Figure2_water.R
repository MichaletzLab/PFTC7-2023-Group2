# --- Soil Moisture GAM ---
gam_iWUE <- mgcv::gam(I(A/gsw) ~ s(Tleaf, by = moist_cat, bs = "cs") + moist_cat,
                      data = dat_environ)

# Prediction grid
pred_grid <- expand.grid(
  Tleaf = seq(min(dat_environ$Tleaf, na.rm = TRUE),
              max(dat_environ$Tleaf, na.rm = TRUE),
              length.out = 500),
  moist_cat = unique(dat_environ$moist_cat)
)

pred <- predict(gam_iWUE, newdata = pred_grid, se.fit = TRUE)

pred_grid <- pred_grid %>%
  mutate(y_pred = pred$fit,
         se = pred$se.fit,
         ymin = y_pred - 1.96 * se,
         ymax = y_pred + 1.96 * se)

# Topt per category
topt_by_bin <- pred_grid %>%
  group_by(moist_cat) %>%
  slice_max(order_by = y_pred, n = 1, with_ties = FALSE)

# Overall GAM (no grouping)
gam_iWUE_all <- mgcv::gam(I(A/gsw) ~ s(Tleaf, bs = "cs"), data = dat_environ)
pred_all <- data.frame(Tleaf = seq(min(dat_environ$Tleaf, na.rm = TRUE),
                                   max(dat_environ$Tleaf, na.rm = TRUE),
                                   length.out = 500))
pred_all$y_pred <- predict(gam_iWUE_all, newdata = pred_all)
overall_topt <- pred_all[which.max(pred_all$y_pred), "Tleaf"]

# Colors
moist_colors <- c("Low" = "grey30",
                  "Medium" = "lightblue3",
                  "High" = "#3333CC")

# Filter for plotting
dat_plot <- dat_environ %>% filter((A/gsw) < 400)
pred_grid_plot <- pred_grid %>% filter(y_pred < 400)

# Plot
iWUE.moist.plot <- ggplot() +
  geom_point(data = dat_plot, aes(x = Tleaf, y = A/gsw, color = moist_cat), alpha = 0.01) +
  geom_ribbon(data = pred_grid_plot,
              aes(x = Tleaf, ymin = ymin, ymax = ymax, fill = moist_cat), alpha = 0.6) +
  geom_line(data = pred_grid_plot, aes(x = Tleaf, y = y_pred, color = moist_cat), size = 1.2) +
  geom_vline(xintercept = overall_topt, linetype = "dashed", size = 1) +
  geom_vline(data = topt_by_bin, aes(xintercept = Tleaf, color = moist_cat),
             linetype = "dotted", size = 1) +
  labs(x = "Leaf temperature (°C)",
       y = expression(iWUE~(mu*mol~CO[2]~mol^{-1}~H[2]*O)),
       color = "Soil moisture",
       fill = "Soil moisture") +
  scale_color_manual(values = moist_colors) +
  scale_fill_manual(values = moist_colors) +
  theme_classic() +
  coord_cartesian(ylim = c(0, 400))
iWUE.moist.plot

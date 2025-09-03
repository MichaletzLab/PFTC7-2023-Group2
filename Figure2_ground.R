# --- 1) Fit GAM with separate smooths for each T2_cat ---
gam_iWUE_T2 <- mgcv::gam(I(A/gsw) ~ s(Tleaf, by = T2_cat, bs = "cs") + T2_cat, 
                         data = dat.T2Temp)

# --- 2) Prediction grid ---
pred_grid_T2 <- expand.grid(
  Tleaf = seq(min(dat.T2Temp$Tleaf, na.rm = TRUE),
              max(dat.T2Temp$Tleaf, na.rm = TRUE),
              length.out = 500),
  T2_cat = levels(dat.T2Temp$T2_cat)
)

pred_T2 <- predict(gam_iWUE_T2, newdata = pred_grid_T2, se.fit = TRUE)

pred_grid_T2 <- pred_grid_T2 %>%
  mutate(
    y_pred = pred_T2$fit,
    se = pred_T2$se.fit,
    ymin = y_pred - 1.96 * se,
    ymax = y_pred + 1.96 * se
  )

# --- Topt per bin ---
topt_by_T2 <- pred_grid_T2 %>%
  group_by(T2_cat) %>%
  slice_max(order_by = y_pred, n = 1, with_ties = FALSE) %>%
  select(T2_cat, Tleaf, y_pred)

# --- Overall GAM (global smooth) ---
gam_T2_all <- mgcv::gam(I(A/gsw) ~ s(Tleaf, bs = "cs"), data = dat.T2Temp)
pred_all_T2 <- data.frame(Tleaf = seq(min(dat.T2Temp$Tleaf, na.rm = TRUE),
                                      max(dat.T2Temp$Tleaf, na.rm = TRUE),
                                      length.out = 500))
pred_all_T2$y_pred <- predict(gam_T2_all, newdata = pred_all_T2)
overall_topt_T2 <- pred_all_T2[which.max(pred_all_T2$y_pred), "Tleaf"]

# --- Define colors ---
T2_colors <- c("Low" = "grey30",
               "Medium" = "lightblue3",
               "High" = "#3333CC")

# --- 3) Filter for display ---
dat_T2_plot <- dat.T2Temp %>% filter((A/gsw) < 400)
pred_grid_T2_plot <- pred_grid_T2 %>% filter(y_pred < 400)

# --- 4) Plot ---
iWUE.Tleaf.T2.plot <- ggplot() +
  geom_point(data = dat_T2_plot,
             aes(x = Tleaf, y = A/gsw, color = T2_cat),
             alpha = 0.01) +
  geom_ribbon(data = pred_grid_T2_plot,
              aes(x = Tleaf, ymin = ymin, ymax = ymax, fill = T2_cat),
              alpha = 0.6) +
  geom_line(data = pred_grid_T2_plot,
            aes(x = Tleaf, y = y_pred, color = T2_cat),
            size = 1.2) +
  geom_vline(xintercept = overall_topt_T2, linetype = "dashed", size = 1) +
  geom_vline(data = topt_by_T2,
             aes(xintercept = Tleaf, color = T2_cat),
             linetype = "dotted", size = 1) +
  labs(
    x = "Leaf temperature (°C)",
    y = expression(iWUE~(mu*mol~CO[2]~mol^{-1}~H[2]*O)),
    color = "Ground temperature (T2)",
    fill = "Ground temperature (T2)"
  ) +
  scale_color_manual(values = T2_colors) +
  scale_fill_manual(values = T2_colors) +
  theme_classic() +
  coord_cartesian(ylim = c(0, 400))

iWUE.Tleaf.T2.plot
